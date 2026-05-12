"""CLS 随机介质输运测试。

参照 CLS_with_VL 分支的 LCLS 测试模式，使用 PyAPITestHarness。
"""

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness

# Path to locally-built OpenMC binary (with stochastic media support)
OPENMC_EXEC = '/home/openmc/build/bin/openmc'


@pytest.fixture
def cls_eigenvalue_model():
    """创建 CLS eigenvalue 测试模型。

    简化球形几何：
    - 颗粒材料: UO2 (5% U235), 密度 10.5 g/cm3
    - 基体材料: Graphite (C0), 密度 1.7 g/cm3
    - 颗粒球面 r=0.02 cm
    - 容器球面 R=1.0 cm, reflective 边界
    - 填充率 phi=0.20
    """
    # Materials
    uo2 = openmc.Material(name='UO2')
    uo2.add_nuclide('U235', 0.05)
    uo2.add_nuclide('U238', 0.95)
    uo2.add_nuclide('O16', 2.0)
    uo2.set_density('g/cm3', 10.5)

    graphite = openmc.Material(name='Graphite')
    graphite.add_nuclide('C0', 1.0)
    graphite.set_density('g/cm3', 1.7)

    # Particle universe (particle cell + remainder cell)
    r = 0.02
    particle_sphere = openmc.Sphere(r=r)
    particle_cell = openmc.Cell(fill=uo2, region=-particle_sphere)
    remainder_cell = openmc.Cell(fill=graphite, region=+particle_sphere)
    particle_universe = openmc.Universe(cells=[particle_cell, remainder_cell])

    # CLS stochastic media
    cls = openmc.CLSMedia(stochastic_media_id=101, name='test_cls')
    cls.pack_fraction = 0.20
    cls.particle_universe = particle_universe
    cls.particle_cell = particle_cell
    cls.matrix_material = graphite

    # Container cell
    R = 1.0
    container_surf = openmc.Sphere(r=R, boundary_type='reflective')
    fuel_cell = openmc.Cell(fill=cls, region=-container_surf)

    root_univ = openmc.Universe(cells=[fuel_cell])
    geometry = openmc.Geometry(root_univ)

    # Settings
    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.batches = 10
    settings.inactive = 5
    settings.particles = 1000
    settings.seed = 1
    source = openmc.IndependentSource(
        space=openmc.stats.Point((0.0, 0.0, 0.0))
    )
    settings.source = source

    materials = openmc.Materials([uo2, graphite])
    return openmc.Model(geometry=geometry, materials=materials,
                        settings=settings)


def test_cls_eigenvalue_smoke(cls_eigenvalue_model, tmp_path):
    """烟囱测试：验证 CLS 模型能完成 eigenvalue 输运不崩溃。"""
    model = cls_eigenvalue_model
    model.export_to_xml(tmp_path)

    openmc.run(cwd=tmp_path, openmc_exec=OPENMC_EXEC)

    # 验证 statepoint 文件存在
    sp_files = list(tmp_path.glob('statepoint.*.h5'))
    assert len(sp_files) == 1, f'Expected 1 statepoint file, found {len(sp_files)}'

    # 通过 OpenMC API 读取 keff
    sp = openmc.StatePoint(sp_files[0])
    assert sp.keff.n > 0, f'keff should be positive, got {sp.keff}'
