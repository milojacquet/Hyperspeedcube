use hypermath::pga::*;
use hypermath::prelude::*;

use super::*;

/// Lua wrapper around a blade in the projective geometric algebra.
#[derive(Debug, Clone)]
pub struct LuaBlade(pub Blade);

impl<'lua> FromLua<'lua> for LuaBlade {
    fn from_lua(value: LuaValue<'lua>, lua: &'lua Lua) -> LuaResult<Self> {
        let hyperbolic = LuaSpace::get(lua)?.hyperbolic();

        if let Ok(m) = cast_userdata(lua, &value) {
            Ok(m)
        } else if let Ok(n) = lua.unpack(value.clone()) {
            Ok(Self(Blade::scalar(LuaNdim::get(lua)?, n, hyperbolic)))
        } else if let Ok(LuaVector(v)) = lua.unpack(value.clone()) {
            let ndim = enforce_ndim(lua, v.ndim())?;
            Ok(Self(Blade::from_vector(ndim, v, hyperbolic)))
        } else if let Ok(LuaHyperplane(h)) = lua.unpack(value.clone()) {
            let ndim = enforce_ndim(lua, h.normal().ndim())?;
            Ok(Self(Blade::from_hyperplane(ndim, &h)))
        } else {
            lua_convert_err(&value, "blade (scalar, vector, point, hyperplane, etc.)")
        }
    }
}

impl LuaUserData for LuaBlade {
    fn add_fields<'lua, F: LuaUserDataFields<'lua, Self>>(fields: &mut F) {
        fields.add_meta_field("type", LuaStaticStr("blade"));

        fields.add_field_method_get("ndim", |_lua, Self(b)| Ok(b.ndim()));
        fields.add_field_method_get("grade", |_lua, Self(b)| Ok(b.grade()));
        fields.add_field_method_get("antigrade", |_lua, Self(b)| Ok(b.antigrade()));

        fields.add_field_method_get("dual", |_lua, Self(this)| Ok(Self(this.dual())));
        fields.add_field_method_get("antidual", |_lua, Self(this)| Ok(Self(this.antidual())));
        fields.add_field_method_get("complement", |_lua, Self(this)| {
            Ok(Self(this.right_complement()))
        });

        fields.add_field_method_get("unit", |_lua, Self(this)| {
            match hypermath::util::try_div(this, this.mag()) {
                Some(v) => Ok(Self(v)),
                None => Err(LuaError::external("cannot normalize zero vector")),
            }
        });
        fields.add_field_method_get("normalize", |_lua, _| {
            Err::<LuaValue<'_>, _>(LuaError::external("use `.unit` instead"))
        });
        fields.add_field_method_get("normalized", |_lua, Self(_)| {
            Err::<LuaValue<'_>, _>(LuaError::external("use `.unit` instead"))
        });

        fields.add_field_method_get("mag2", |_lua, Self(this)| Ok(this.mag2()));
        fields.add_field_method_get("mag", |_lua, Self(this)| Ok(this.mag()));

        fields.add_field_method_get("bulk", |_lua, Self(this)| Ok(Self(this.bulk())));
        fields.add_field_method_get("weight", |_lua, Self(this)| Ok(Self(this.weight())));
    }

    fn add_methods<'lua, M: LuaUserDataMethods<'lua, Self>>(methods: &mut M) {
        methods.add_meta_method(LuaMetaMethod::ToString, |_lua, Self(b), ()| {
            if let Some(p) = b.to_point() {
                Ok(format!("point{p}"))
            } else if let Some(v) = b.to_vector() {
                Ok(format!("vec{v}"))
            } else {
                Ok(format!("blade({b})"))
            }
        });

        methods.add_meta_method("cross", |_lua, Self(this), Self(other)| {
            Ok(Blade::cross_product_3d(&this, &other).map(Self))
        });

        methods.add_method("projected_to", |_lua, Self(this), Self(other)| {
            Ok(this.orthogonal_projection_to(&other).map(Self))
        });
        methods.add_method("rejected_from", |_lua, Self(this), Self(other)| {
            Ok(this.orthogonal_rejection_from(&other).map(Self))
        });

        methods.add_method("wedge", |_lua, Self(this), Self(other)| {
            Ok(Blade::wedge(this, &other).map(Self))
        });
        methods.add_method("antiwedge", |_lua, Self(this), Self(other)| {
            Ok(Blade::antiwedge(this, &other).map(Self))
        });
        methods.add_method("dot", |_lua, Self(this), Self(other)| {
            Ok(Blade::dot(this, &other))
        });
        methods.add_method("antidot", |_lua, Self(this), Self(other)| {
            Ok(Blade::antidot(this, &other).map(|term| Self(Blade::from_term(this.ndim(), term))))
        });

        methods.add_method(
            "cross",
            |_lua, Self(this), Self(other)| match Blade::cross_product_3d(&this, &other) {
                Some(result) => Ok(LuaBlade(result)),
                None => Err(LuaError::external(
                    "cross product is only allowed between vectors in 3D",
                )),
            },
        );

        // blade + blade
        methods.add_meta_function(LuaMetaMethod::Add, |_lua, (Self(lhs), Self(rhs))| {
            Ok(Self(lhs + rhs))
        });
        // blade - blade
        methods.add_meta_function(LuaMetaMethod::Sub, |_lua, (Self(lhs), Self(rhs))| {
            Ok(Self(lhs - rhs))
        });

        // blade * number; number * blade
        methods.add_meta_function(LuaMetaMethod::Mul, |lua, args: LuaMultiValue<'_>| {
            if let Ok((Self(v), a)) = lua.unpack_multi(args.clone()) {
                let a: Float = a;
                Ok(Self(v * a))
            } else if let Ok((a, Self(v))) = lua.unpack_multi(args.clone()) {
                let a: Float = a;
                Ok(Self(v * a))
            } else {
                let [a, b]: [LuaValue<'_>; 2] = lua.unpack_multi(args)?;
                Err(LuaError::external(format!(
                    "cannot multiply {} and {}",
                    a.type_name(),
                    b.type_name(),
                )))
            }
        });

        // blade / number
        methods.add_meta_method(LuaMetaMethod::Div, |_lua, Self(lhs), rhs| {
            Ok(hypermath::util::try_div(lhs, rhs).map(Self))
        });

        // -blade
        methods.add_meta_method(LuaMetaMethod::Unm, |_lua, Self(b), ()| Ok(Self(-b)));

        // blade ^ blade
        methods.add_meta_function(LuaMetaMethod::Pow, |_lua, (Self(a), Self(b))| {
            Ok(Blade::wedge(&a, &b).map(Self))
        });
        // blade & blade
        methods.add_meta_function(LuaMetaMethod::BAnd, |_lua, (Self(a), Self(b))| {
            Ok(Blade::antiwedge(&a, &b).map(Self))
        });
        // ~blade
        methods.add_meta_function(LuaMetaMethod::BNot, |_lua, Self(this)| {
            Ok(Self(this.right_complement()))
        });

        // blade == blade
        methods.add_meta_function(LuaMetaMethod::Eq, |_lua, (Self(a), Self(b))| {
            Ok(approx_eq(&a, &b))
        });

        // blade[index]
        methods.add_meta_method(
            LuaMetaMethod::Index,
            |lua, Self(this), arg: LuaValue<'_>| match lua.unpack(arg) {
                Ok(LuaMultivectorIndex { axes, sign, .. }) => Ok(this.get(axes).map(|&x| x * sign)),
                Err(_) => Ok(None),
            },
        );

        // We do not support `LuaMetaMethod::NewIndex` because this can be used
        // to mutate aliased blades, which is very confusing.
        methods.add_meta_method(
            LuaMetaMethod::NewIndex,
            |_lua, Self(_), _: LuaMultiValue<'_>| -> LuaResult<()> {
                Err(LuaError::external(
                    "mutation of blades is not allowed. \
                     construct a new blade instead.",
                ))
            },
        );

        // pairs(Vector)
        methods.add_meta_function(LuaMetaMethod::Pairs, |lua, this: LuaValue<'_>| {
            match lua.unpack(this.clone()) {
                Ok(LuaVector(_)) => {
                    let vector_iter =
                        lua.create_function(|_lua, (LuaVector(v), LuaVectorIndex(i))| {
                            if i < v.ndim() {
                                Ok((Some(LuaVectorIndex(i + 1)), Some(v[i])))
                            } else {
                                Ok((None, None))
                            }
                        })?;

                    Ok((vector_iter, this, LuaVectorIndex(0)))
                }

                Err(_) => Err(LuaError::external(
                    "iteration is only supported for blades representing a point or vector",
                )),
            }
        });
    }
}

impl LuaBlade {
    /// Constructs a blade representing a vector.
    pub fn from_vector(lua: &Lua, v: impl VectorRef) -> LuaResult<Self> {
        let hyperbolic = LuaSpace::get(lua)?.hyperbolic();
        Ok(Self(Blade::from_vector(LuaNdim::get(lua)?, v, hyperbolic)))
    }
    /// Constructs a blade representing a point.
    pub fn from_point(lua: &Lua, v: impl VectorRef) -> LuaResult<Self> {
        let hyperbolic = LuaSpace::get(lua)?.hyperbolic();
        Ok(Self(Blade::from_point(LuaNdim::get(lua)?, v, hyperbolic)))
    }
}

impl TransformByMotor for LuaBlade {
    fn transform_by(&self, m: &Motor) -> Self {
        Self(self.0.transform_by(m))
    }
}

fn enforce_ndim(lua: &Lua, ndim: u8) -> LuaResult<u8> {
    let expected_ndim = LuaNdim::get(lua)?;
    if ndim <= expected_ndim {
        Ok(ndim)
    } else {
        Err(LuaError::external(format!(
            "cannot construct {ndim}D object in {expected_ndim}D space",
        )))
    }
}
