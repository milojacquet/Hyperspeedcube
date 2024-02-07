use super::*;

pipeline!(pub(in crate::gfx) struct Pipeline {
    type = wgpu::RenderPipeline;

    struct Bindings<'a> {
        sticker_colors_texture: &'a wgpu::TextureView = pub(FRAGMENT) bindings::STICKER_COLORS_TEXTURE,
        special_colors_texture: &'a wgpu::TextureView = pub(FRAGMENT) bindings::SPECIAL_COLORS_TEXTURE,
        ids_texture:            &'a wgpu::TextureView = pub(FRAGMENT) bindings::IDS_TEXTURE,
        polygon_color_ids:      &'a wgpu::Buffer = pub(FRAGMENT) bindings::POLYGON_COLOR_IDS,
        edges:                  &'a wgpu::Buffer = pub(FRAGMENT) bindings::EDGE_VERTS,
        vertex_3d_positions:    &'a wgpu::Buffer = pub(FRAGMENT) bindings::VERTEX_3D_POSITIONS_READONLY,
        lighting_params:        &'a wgpu::Buffer = pub(FRAGMENT) bindings::LIGHTING_PARAMS,
        composite_params:       &'a wgpu::Buffer = pub(FRAGMENT) bindings::COMPOSITE_PARAMS,

        projection_params:      &'a wgpu::Buffer = pub(FRAGMENT) bindings::PROJECTION_PARAMS,
        view_params:            &'a wgpu::Buffer = pub(FRAGMENT) bindings::VIEW_PARAMS,
        target_size:            &'a wgpu::Buffer = pub(FRAGMENT) bindings::TARGET_SIZE,
    }

    let pipeline_descriptor = RenderPipelineDescriptor {
        label: "render_composite_puzzle",
        vertex_entry_point: "uv_vertex",
        fragment_entry_point: "render_composite_puzzle_fragment",
        vertex_buffers: &[UvVertex::LAYOUT],
        primitive: wgpu::PrimitiveState {
            topology: wgpu::PrimitiveTopology::TriangleStrip,
            ..Default::default()
        },
        fragment_target: Some(wgpu::ColorTargetState {
            format: wgpu::TextureFormat::Rgba8UnormSrgb,
            blend: Some(wgpu::BlendState {
                color: blend_component!(Add(src * Constant, dst * One)),
                alpha: blend_component!(Add(src * Constant, dst * One)),
            }),
            write_mask: wgpu::ColorWrites::ALL,
        }),
        ..Default::default()
    };
});

pub(in crate::gfx) struct PassParams<'tex> {
    pub clear: bool,
    pub target: &'tex wgpu::TextureView,
}
impl<'pass> PassParams<'pass> {
    pub fn begin_pass(self, encoder: &'pass mut wgpu::CommandEncoder) -> wgpu::RenderPass<'pass> {
        encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("render_composite_puzzle"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view: self.target,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: match self.clear {
                        true => wgpu::LoadOp::Clear(wgpu::Color::TRANSPARENT),
                        false => wgpu::LoadOp::Load,
                    },
                    store: wgpu::StoreOp::Store,
                },
            })],
            ..Default::default()
        })
    }
}
