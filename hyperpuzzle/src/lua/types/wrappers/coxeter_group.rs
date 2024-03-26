use hypershape::CoxeterGroup;
use std::str::FromStr;

use super::*;

/// Conversion wrapper for a Coxeter group.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct LuaCoxeterGroupSymbol(pub CoxeterGroup);
impl<'lua> FromLua<'lua> for LuaCoxeterGroupSymbol {
    fn from_lua(lua_value: LuaValue<'lua>, lua: LuaContext<'lua>) -> LuaResult<Self> {
        lua_convert!(match (lua, &lua_value, "coxeter group") {
            LuaValue::String(s) => s
                .to_str()?
                .parse()
                .map(|LuaCoxeterGroupSymbol(c)| LuaCoxeterGroupSymbol(c)),
        })
    }
}

fn split_alpha_number(s: &str) -> Option<(&str, u8)> {
    let (num_index, _) = s.match_indices(char::is_numeric).next()?;
    let num = s[num_index..].parse().ok()?;
    Some((&s[0..num_index], num))
}

impl FromStr for LuaCoxeterGroupSymbol {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.to_lowercase();
        match &s[..] {
            "e6" => Ok(LuaCoxeterGroupSymbol(CoxeterGroup::E6)),
            "e7" => Ok(LuaCoxeterGroupSymbol(CoxeterGroup::E7)),
            "e8" => Ok(LuaCoxeterGroupSymbol(CoxeterGroup::E8)),
            "f4" => Ok(LuaCoxeterGroupSymbol(CoxeterGroup::F4)),
            "g2" => Ok(LuaCoxeterGroupSymbol(CoxeterGroup::G2)),
            "h2" => Ok(LuaCoxeterGroupSymbol(CoxeterGroup::H2)),
            "h3" => Ok(LuaCoxeterGroupSymbol(CoxeterGroup::H3)),
            "h4" => Ok(LuaCoxeterGroupSymbol(CoxeterGroup::H4)),
            _ => match split_alpha_number(&s) {
                Some(("a", n)) => Ok(LuaCoxeterGroupSymbol(CoxeterGroup::A(n))),
                Some(("b" | "c" | "bc", n)) => Ok(LuaCoxeterGroupSymbol(CoxeterGroup::B(n))),
                Some(("d", n)) => Ok(LuaCoxeterGroupSymbol(CoxeterGroup::D(n))),
                Some(("i", n)) => Ok(LuaCoxeterGroupSymbol(CoxeterGroup::I(n))),
                _ => Err("invalid coxeter group".to_owned()),
            },
        }
    }
}
