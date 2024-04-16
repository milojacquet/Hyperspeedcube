use itertools::Itertools;
use parking_lot::{Mutex, MutexGuard};
use std::{borrow::Cow, collections::HashMap, hash::Hash, sync::Arc};

use hypermath::prelude::*;

use crate::builder::{CustomOrdering, NamingScheme};

use super::*;

#[derive(Debug)]
pub struct LuaDbEntry<I, D> {
    pub id: I,
    pub db: Arc<Mutex<D>>,
}
impl<I: Clone, D> Clone for LuaDbEntry<I, D> {
    fn clone(&self) -> Self {
        Self {
            id: self.id.clone(),
            db: Arc::clone(&self.db),
        }
    }
}

impl<'lua, I, D> FromLua<'lua> for LuaDbEntry<I, D>
where
    Self: 'static + LuaUserData + Clone,
{
    fn from_lua(value: LuaValue<'lua>, lua: &'lua Lua) -> LuaResult<Self> {
        cast_userdata(lua, &value)
    }
}

/// Database of Lua values referenced using some sort of unique ID.
pub trait LuaIdDatabase<'lua, I>: 'static + Sized
where
    I: 'static + Clone,
    LuaDbEntry<I, Self>: LuaUserData,
{
    const ELEMENT_NAME_SINGULAR: &'static str;
    const ELEMENT_NAME_PLURAL: &'static str;

    /// Converts the ID of an entry to a [`LuaValue`].
    fn wrap_id(&self, id: I) -> LuaDbEntry<I, Self> {
        let db = self.db_arc();
        LuaDbEntry { id, db }
    }
    /// Converts a [`LuaValue`] to an entry ID, or returns an error if no such
    /// entry exists.
    fn value_to_id(&self, lua: &'lua Lua, value: LuaValue<'lua>) -> LuaResult<I>;

    /// Converts a [`LuaValue`] to an entry ID if it is a [`LuaDbEntry`]
    /// userdata value, or returns `None` if it is not.
    fn value_to_id_by_userdata(
        &self,
        lua: &'lua Lua,
        value: &LuaValue<'lua>,
    ) -> Option<LuaResult<I>>
    where
        LuaDbEntry<I, Self>: LuaUserData,
    {
        let LuaDbEntry { id, db } = cast_userdata(lua, value).ok()?;
        Some(match Arc::ptr_eq(&self.db_arc(), &db) {
            true => Ok(id),
            false => Err(LuaError::external(
                "cannot operate on entries from a different database",
            )),
        })
    }

    /// Returns an `Arc` reference to the DB.
    fn db_arc(&self) -> Arc<Mutex<Self>>;
    /// Returns the number of entries in the database.
    fn db_len(&self) -> usize;
    /// Returns a list of IDs in the database, ideally in some canonical order.
    fn ids_in_order(&self) -> Cow<'_, [I]>;

    fn mapping_from_value<T: FromLua<'lua>>(
        &self,
        lua: &'lua Lua,
        mapping_value: LuaValue<'lua>,
    ) -> LuaResult<Vec<(I, T)>> {
        match mapping_value {
            LuaValue::Table(t) => t
                .pairs()
                .map(|pair| {
                    let (id, new_value) = pair?;
                    Ok((self.value_to_id(lua, id)?, new_value))
                })
                .try_collect(),

            LuaValue::Function(f) => self
                .ids_in_order()
                .iter()
                .map(|id| {
                    let new_value = f.call(self.wrap_id(id.clone()))?;
                    Ok((id.clone(), new_value))
                })
                .try_collect(),

            _ => return lua_convert_err(&mapping_value, "table or function"),
        }
    }

    /// Defines the following methods:
    /// - `__index` (metamethod)
    /// - `__len` (metamethod)
    fn add_db_metamethods<T: 'static, M: LuaUserDataMethods<'lua, T>>(
        methods: &mut M,
        lock: fn(&T) -> MutexGuard<'_, Self>,
    ) {
        methods.add_meta_method(LuaMetaMethod::Index, move |lua, this, index| {
            let db = lock(this);
            match db.value_to_id(lua, index) {
                Ok(id) => Ok(Some(db.wrap_id(id))),
                Err(_) => Ok(None),
            }
        });
        methods.add_meta_method(LuaMetaMethod::Len, move |_lua, this, ()| {
            let db = lock(this);
            Ok(db.db_len())
        });
    }
}

/// Extension of [`LuaIdDatabase`] to support naming elements.
pub trait LuaNamedIdDatabase<'lua, I>: LuaIdDatabase<'lua, I>
where
    I: 'static + Clone + Hash + Eq,
    LuaDbEntry<I, Self>: LuaUserData,
{
    /// Returns a reference to the naming scheme of the database.
    fn names(&self) -> &NamingScheme<I>;
    /// Returns a mutable reference to the naming scheme of the database.
    fn names_mut(&mut self) -> &mut NamingScheme<I>;

    /// Converts a [`LuaValue`] to an entry ID if it is a string containing an
    /// element name, or returns `None` if it is not.
    #[must_use]
    fn value_to_id_by_name(&self, _lua: &'lua Lua, value: &LuaValue<'lua>) -> Option<LuaResult<I>> {
        let s = value.as_str()?;
        Some(match self.names().names_to_ids().get(s) {
            Some(id) => Ok(id.clone()),
            None => Err(LuaError::external(format!("no entry named {s:?}"))),
        })
    }

    /// Renames all elements according to a table or function.
    fn rename_all(
        &mut self,
        lua: &'lua Lua,
        new_names_value: LuaValue<'lua>,
    ) -> LuaResult<NamingScheme<I>> {
        // We need to rename all the entries at once, so just construct a new
        // naming scheme from scratch.
        let mut new_names = NamingScheme::new();

        // First, assemble a list of all the renames that need to happen.
        let kv_pairs: Vec<(I, Option<String>)> = self.mapping_from_value(lua, new_names_value)?;

        // Set the new names.
        for (k, v) in kv_pairs {
            new_names.set(k, v).into_lua_err()?;
        }

        Ok(new_names)
    }

    /// Defines the following methods on a database:
    /// - `rename`
    fn add_named_db_methods<T: 'static, M: LuaUserDataMethods<'lua, T>>(
        methods: &mut M,
        lock: fn(&T) -> MutexGuard<'_, Self>,
    ) {
        methods.add_method("rename", move |lua, this, new_names| {
            let mut db = lock(this);
            *db.names_mut() = db.rename_all(lua, new_names)?;
            Ok(())
        });
    }

    /// Defines the following fields on a database entry:
    /// - `name`
    fn add_named_db_entry_fields<F: LuaUserDataFields<'lua, LuaDbEntry<I, Self>>>(fields: &mut F) {
        fields.add_field_method_get("name", |_lua, this| {
            let db = this.db.lock();
            Ok(db.names().get(this.id.clone()))
        });
        fields.add_field_method_set("name", |_lua, this, new_name| {
            let mut db = this.db.lock();
            db.names_mut().set(this.id.clone(), new_name).into_lua_err()
        });
    }
}

/// Extension of [`LuaIdDatabase`] to enforce a total ordering on entries. Also,
/// the ID must be an [`IndexNewtype`].
pub trait LuaOrderedIdDatabase<'lua, I>: LuaIdDatabase<'lua, I>
where
    I: 'static + IndexNewtype,
    LuaDbEntry<I, Self>: LuaUserData,
{
    /// Returns a reference to the custom ordering of entries in the database.
    fn ordering(&self) -> &CustomOrdering<I>;
    /// Returns a mutable reference to the custom ordering of entries in the
    /// database.
    fn ordering_mut(&mut self) -> &mut CustomOrdering<I>;

    /// Converts a [`LuaValue`] to an entry ID if it is an index, or returns
    /// `None` if it is not.
    fn value_to_id_by_index(&self, lua: &'lua Lua, value: &LuaValue<'lua>) -> Option<LuaResult<I>> {
        let LuaIndex(i) = lua.unpack(value.clone()).ok()?;
        Some(match self.ordering().ids_in_order().get(i) {
            Some(&id) => Ok(id),
            None => Err(LuaError::external(if self.db_len() == 1 {
                format!(
                    "index {} is out of range; there is only 1 {}",
                    i + 1,
                    Self::ELEMENT_NAME_SINGULAR,
                )
            } else {
                format!(
                    "index {} is out of range; there are only {} {}",
                    i + 1,
                    self.db_len(),
                    Self::ELEMENT_NAME_PLURAL,
                )
            })),
        })
    }

    /// Reorders all elements according to a table or function.
    fn reorder_all(&mut self, lua: &'lua Lua, new_ordering_value: LuaValue<'lua>) -> LuaResult<()> {
        let new_order_keys: HashMap<I, f64> = if let LuaValue::Table(t) = new_ordering_value {
            t.sequence_values()
                .enumerate()
                .map(|(i, elem)| LuaResult::Ok((self.value_to_id(lua, elem?)?, i as f64)))
                .try_collect()?
        } else {
            self.mapping_from_value(lua, new_ordering_value)?
                .into_iter()
                .collect()
        };

        // By default, leave unspecified entries in the same order at the end.
        // This sort is guaranteed to be stable.
        let mut new_ordering: Vec<I> = self.ordering().ids_in_order().to_vec();
        new_ordering.sort_by(|a, b| {
            f64::total_cmp(
                new_order_keys.get(a).unwrap_or(&f64::INFINITY),
                new_order_keys.get(b).unwrap_or(&f64::INFINITY),
            )
        });

        // We will apply the new ordering all at once.
        let current_ordering = self.ordering_mut();
        for (index, id) in new_ordering.into_iter().enumerate() {
            current_ordering.swap_to_index(id, index).into_lua_err()?;
        }

        Ok(())
    }

    fn swap(&mut self, lua: &'lua Lua, i: LuaValue<'lua>, j: LuaValue<'lua>) -> LuaResult<()> {
        let i = self.value_to_id(lua, i)?;
        let j = self.value_to_id(lua, j)?;
        self.ordering_mut().swap(i, j);
        Ok(())
    }

    /// Defines the following methods on a database:
    /// - `swap`
    /// - `reorder`
    fn add_ordered_db_methods<T: 'static, M: LuaUserDataMethods<'lua, T>>(
        methods: &mut M,
        lock: fn(&T) -> MutexGuard<'_, Self>,
    ) {
        methods.add_method("swap", move |lua, this, (i, j)| {
            let mut db = lock(this);
            db.swap(lua, i, j)
        });

        methods.add_method("reorder", move |lua, this, new_ordering| {
            let mut db = lock(this);
            db.reorder_all(lua, new_ordering)
        });
    }

    /// Defines the following fields on a database entry:
    /// - `index`
    fn add_ordered_db_entry_fields<F: LuaUserDataFields<'lua, LuaDbEntry<I, Self>>>(
        fields: &mut F,
    ) {
        fields.add_field_method_get("index", |_lua, this| {
            let db = this.db.lock();
            db.ordering().get_index(this.id).into_lua_err()
        });
        fields.add_field_method_set("index", |lua, this, new_index| {
            let mut db = this.db.lock();
            let new_index = db.value_to_id(lua, new_index)?;
            db.ordering_mut().shift_to(this.id, new_index);
            Ok(())
        });
    }
}