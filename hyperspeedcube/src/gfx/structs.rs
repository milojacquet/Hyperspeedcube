//! Structs shared between the CPU and GPU (vertices, uniforms, etc.).

#[repr(C)]
#[derive(Debug, Copy, Clone, PartialEq, bytemuck::NoUninit, bytemuck::Zeroable)]
pub(super) struct GfxLightingParams {
    pub dir: [f32; 3],
    pub ambient: f32,
}

#[repr(C)]
#[derive(Debug, Copy, Clone, PartialEq, bytemuck::NoUninit, bytemuck::Zeroable)]
pub(super) struct GfxProjectionParams {
    pub facet_shrink: f32,
    pub sticker_shrink: f32,
    pub piece_explode: f32,

    pub w_factor_4d: f32,
    pub w_factor_3d: f32,
    pub fov_signum: f32,
}

#[repr(C)]
#[derive(Debug, Default, Copy, Clone, PartialEq, bytemuck::NoUninit, bytemuck::Zeroable)]
pub(super) struct GfxCompositeParams {
    pub outline_radius: u32,
}

#[repr(C)]
#[derive(Debug, Default, Copy, Clone, PartialEq, bytemuck::NoUninit, bytemuck::Zeroable)]
pub(super) struct UvVertex {
    pub position: [f32; 2],
    pub uv: [f32; 2],
}
impl UvVertex {
    const fn new(position: [f32; 2], uv: [f32; 2]) -> Self {
        Self { position, uv }
    }
}
impl UvVertex {
    pub const LAYOUT: wgpu::VertexBufferLayout<'static> = wgpu::VertexBufferLayout {
        array_stride: std::mem::size_of::<Self>() as wgpu::BufferAddress,
        step_mode: wgpu::VertexStepMode::Vertex,
        attributes: &wgpu::vertex_attr_array![
            0 => Float32x2,
            1 => Float32x2,
        ],
    };
    pub const SQUARE: [Self; 4] = [
        UvVertex::new([-1.0, 1.0], [0.0, 0.0]),
        UvVertex::new([1.0, 1.0], [1.0, 0.0]),
        UvVertex::new([-1.0, -1.0], [0.0, 1.0]),
        UvVertex::new([1.0, -1.0], [1.0, 1.0]),
    ];
}

#[repr(C)]
#[derive(Debug, Copy, Clone, PartialEq, bytemuck::NoUninit, bytemuck::Zeroable)]
pub struct GfxViewParams {
    pub scale: [f32; 2],
    pub align: [f32; 2],

    pub clip_4d_backfaces: i32,     // bool
    pub clip_4d_behind_camera: i32, // bool
}