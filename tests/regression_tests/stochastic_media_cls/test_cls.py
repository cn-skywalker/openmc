"""CLS 随机介质基础集成测试。"""

import os
import tempfile
import openmc
import pytest


def make_particle_universe(uo2, graphite, r=0.02):
    """创建标准 particle universe（颗粒 cell + 剩余 cell）。"""
    particle_sphere = openmc.Sphere(r=r)
    particle_cell = openmc.Cell(fill=uo2, region=-particle_sphere)
    remainder_cell = openmc.Cell(fill=graphite, region=+particle_sphere)
    return openmc.Universe(cells=[particle_cell, remainder_cell]), particle_cell


def make_materials():
    """创建测试用材料。"""
    uo2 = openmc.Material(name='UO2')
    uo2.add_nuclide('U235', 0.05)
    uo2.add_nuclide('U238', 0.95)
    uo2.add_nuclide('O16', 2.0)
    uo2.set_density('g/cm3', 10.5)

    graphite = openmc.Material(name='Graphite')
    graphite.add_nuclide('C0', 1.0)
    graphite.set_density('g/cm3', 1.7)

    return uo2, graphite


def make_settings():
    """创建测试用 settings。"""
    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.batches = 10
    settings.particles = 1000
    src = openmc.IndependentSource()
    src.space = openmc.stats.Point((0.0, 0.0, 0.0))
    src.energy = openmc.stats.Discrete([14.0e6], [1.0])
    src.angle = openmc.stats.Isotropic()
    settings.source = src
    return settings


def test_cls_basic():
    """测试 CLS 随机介质基本功能：创建几何、导出 XML、无崩溃。"""
    uo2, graphite = make_materials()
    particle_univ, particle_cell = make_particle_universe(uo2, graphite)

    cls = openmc.CLSMedia(stochastic_media_id=101, name='test_cls')
    cls.pack_fraction = 0.20
    cls.particle_universe = particle_univ
    cls.particle_cell = particle_cell
    cls.matrix_material = graphite

    R = 1.0
    container_surf = openmc.Sphere(r=R, boundary_type='vacuum')
    fuel_cell = openmc.Cell(fill=cls, region=-container_surf)

    root_univ = openmc.Universe(cells=[fuel_cell])
    geometry = openmc.Geometry(root_univ)

    settings = make_settings()
    materials = openmc.Materials([uo2, graphite])
    model = openmc.Model(geometry=geometry, materials=materials, settings=settings)

    with tempfile.TemporaryDirectory() as tmpdir:
        model.export_to_xml(tmpdir)
        assert os.path.isfile(os.path.join(tmpdir, 'geometry.xml'))
        assert os.path.isfile(os.path.join(tmpdir, 'materials.xml'))
        assert os.path.isfile(os.path.join(tmpdir, 'settings.xml'))


def test_cls_pack_fraction_validation():
    """测试 CLSMedia pack_fraction 参数验证。"""
    cls = openmc.CLSMedia(stochastic_media_id=2)

    cls.pack_fraction = 0.3
    assert cls.pack_fraction == 0.3

    cls.pack_fraction = 0.0
    assert cls.pack_fraction == 0.0

    cls.pack_fraction = 1.0
    assert cls.pack_fraction == 1.0

    with pytest.raises(ValueError):
        cls.pack_fraction = -0.1

    with pytest.raises(ValueError):
        cls.pack_fraction = 1.5


def test_cls_xml_roundtrip():
    """测试 CLSMedia XML 序列化。"""
    uo2 = openmc.Material(name='UO2')
    uo2.add_nuclide('U235', 1.0)
    uo2.set_density('g/cm3', 10.5)

    graphite = openmc.Material(name='Graphite')
    graphite.add_nuclide('C0', 1.0)
    graphite.set_density('g/cm3', 1.7)

    particle_univ, particle_cell = make_particle_universe(uo2, graphite, r=0.02)

    cls = openmc.CLSMedia(stochastic_media_id=5, name='test_xml')
    cls.pack_fraction = 0.4
    cls.particle_universe = particle_univ
    cls.particle_cell = particle_cell
    cls.matrix_material = graphite

    elem = cls.to_xml_element()
    assert elem.tag == 'cls_media'
    assert elem.get('id') == '5'
    assert elem.get('name') == 'test_xml'
    assert elem.get('pack_fraction') == '0.4'
    assert elem.get('particle_universe') == str(particle_univ.id)
    assert elem.get('particle_cell') == str(particle_cell.id)
    assert elem.get('matrix_material') == str(graphite.id)


def test_cls_cell_fill_type():
    """测试 Cell 填充 CLSMedia 后的 fill_type。"""
    cls = openmc.CLSMedia(stochastic_media_id=10)
    cell = openmc.Cell(fill=cls)
    assert cell.fill_type == 'stochastic'
    assert cell.fill is cls


def test_cls_seed_property():
    """测试 CLSMedia seed 属性的设置和验证。"""
    cls = openmc.CLSMedia(stochastic_media_id=20)

    # 默认值 None
    assert cls.seed is None

    # 设置正值
    cls.seed = 42
    assert cls.seed == 42

    # 设置零值（允许）
    cls.seed = 0
    assert cls.seed == 0

    # 类型错误
    with pytest.raises(TypeError):
        cls.seed = 3.14

    # 负数错误
    with pytest.raises(ValueError):
        cls.seed = -1


def test_cls_seed_xml_output():
    """测试 CLSMedia seed 属性的 XML 输出。"""
    uo2, graphite = make_materials()
    particle_univ, particle_cell = make_particle_universe(uo2, graphite)

    # 有 seed
    cls = openmc.CLSMedia(stochastic_media_id=30)
    cls.pack_fraction = 0.3
    cls.particle_universe = particle_univ
    cls.particle_cell = particle_cell
    cls.matrix_material = graphite
    cls.seed = 42

    elem = cls.to_xml_element()
    assert elem.get('seed') == '42'

    # 无 seed
    cls2 = openmc.CLSMedia(stochastic_media_id=31)
    cls2.pack_fraction = 0.3
    cls2.particle_universe = particle_univ
    cls2.particle_cell = particle_cell
    cls2.matrix_material = graphite

    elem2 = cls2.to_xml_element()
    assert elem2.get('seed') is None
