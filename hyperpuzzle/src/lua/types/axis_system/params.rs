use std::{
    collections::{BTreeMap, HashMap},
    sync::Arc,
};

use hypershape::Space;
use parking_lot::Mutex;

use crate::{
    builder::AxisSystemBuilder,
    library::{Cached, LibraryDb, LibraryFile, LibraryFileLoadResult, LibraryObjectParams},
};

use super::*;

#[derive(Debug)]
pub struct AxisSystemParams {
    pub id: Option<String>,

    pub ndim: LuaNdim,
    pub symmetry: Option<LuaSymmetry>,

    user_build_fn: LuaRegistryKey,
}

impl<'lua> FromLua<'lua> for AxisSystemParams {
    fn from_lua(value: LuaValue<'lua>, lua: &'lua Lua) -> LuaResult<Self> {
        let table = lua.unpack(value)?;

        let ndim: LuaNdim;
        let symmetry: Option<LuaSymmetry>;
        let build: LuaFunction<'_>;
        unpack_table!(lua.unpack(table {
            ndim,
            symmetry,
            build,
        }));

        Ok(AxisSystemParams {
            id: None,

            ndim,
            symmetry,

            user_build_fn: lua.create_registry_value(build)?,
        })
    }
}

impl LibraryObjectParams for AxisSystemParams {
    const NAME: &'static str = "axis system";

    type Constructed = Arc<Mutex<AxisSystemBuilder>>;

    fn get_file_map(lib: &LibraryDb) -> &BTreeMap<String, Arc<LibraryFile>> {
        &lib.axis_systems
    }
    fn get_id_map_within_file(
        result: &mut LibraryFileLoadResult,
    ) -> &mut HashMap<String, Cached<Self>> {
        &mut result.axis_systems
    }

    fn new_constructed(space: &Arc<Mutex<Space>>) -> LuaResult<Self::Constructed> {
        Ok(AxisSystemBuilder::new(None, Arc::clone(space)))
    }
    fn clone_constructed(
        existing: &Self::Constructed,
        space: &Arc<Mutex<Space>>,
    ) -> LuaResult<Self::Constructed> {
        existing.lock().clone(space).into_lua_err()
    }
    fn build(&self, lua: &Lua, space: &Arc<Mutex<Space>>) -> LuaResult<Self::Constructed> {
        let LuaNdim(self_ndim) = self.ndim;
        let space_ndim = space.lock().ndim();
        if space_ndim != self_ndim {
            return Err(LuaError::external(format!(
                "axis system requires {self_ndim}D space but was given {space_ndim}D space",
            )));
        }

        let axis_builder = AxisSystemBuilder::new(self.id.clone(), Arc::clone(space));

        axis_builder.lock().symmetry = self.symmetry.clone().map(|sym| sym.schlafli);

        let () = LuaSpace(Arc::clone(space)).with_this_as_global_space(lua, || {
            lua.registry_value::<LuaFunction<'_>>(&self.user_build_fn)?
                .call(LuaAxisSystem(Arc::clone(&axis_builder)))
                .context("error executing axis definition")
        })?;

        Ok(axis_builder)
    }
}