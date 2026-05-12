"""CLS 随机介质 verbose 事件验证测试。

使用 settings.verbosity=10 和 settings.trace 追踪粒子，
验证中子事件序列符合 CLS 物理预期流程。
"""

import subprocess
import openmc

OPENMC_EXEC = '/home/openmc/build/bin/openmc'


def make_cls_verbose_model():
    """创建用于 verbose 事件验证的 CLS 模型。"""
    uo2 = openmc.Material(name='UO2')
    uo2.add_nuclide('U235', 0.05)
    uo2.add_nuclide('U238', 0.95)
    uo2.add_nuclide('O16', 2.0)
    uo2.set_density('g/cm3', 10.5)

    graphite = openmc.Material(name='Graphite')
    graphite.add_nuclide('C0', 1.0)
    graphite.set_density('g/cm3', 1.7)

    r = 0.02
    particle_sphere = openmc.Sphere(r=r)
    particle_cell = openmc.Cell(fill=uo2, region=-particle_sphere)
    remainder_cell = openmc.Cell(fill=graphite, region=+particle_sphere)
    particle_universe = openmc.Universe(cells=[particle_cell, remainder_cell])

    cls = openmc.CLSMedia(stochastic_media_id=101, name='test_cls')
    cls.pack_fraction = 0.20
    cls.particle_universe = particle_universe
    cls.particle_cell = particle_cell
    cls.matrix_material = graphite

    R = 1.0
    container_surf = openmc.Sphere(r=R, boundary_type='reflective')
    fuel_cell = openmc.Cell(fill=cls, region=-container_surf)
    root_univ = openmc.Universe(cells=[fuel_cell])
    geometry = openmc.Geometry(root_univ)

    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.batches = 10
    settings.inactive = 5
    settings.particles = 50
    settings.seed = 1
    settings.verbosity = 10
    settings.trace = (1, 1, 1)
    source = openmc.IndependentSource(
        space=openmc.stats.Point((0.0, 0.0, 0.0)))
    settings.source = source

    materials = openmc.Materials([uo2, graphite])
    return openmc.Model(geometry=geometry, materials=materials, settings=settings)


def test_cls_verbose_physics_flow(tmp_path):
    """通过 verbose 输出验证 CLS 中子完整物理流程。"""
    model = make_cls_verbose_model()
    model.export_to_xml(tmp_path)

    import os
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '1'

    result = subprocess.run(
        [OPENMC_EXEC], cwd=tmp_path,
        capture_output=True, text=True, env=env)
    output = result.stdout + result.stderr

    # 阶段 1: 基础完整性
    assert result.returncode == 0, (
        f"OpenMC crashed (rc={result.returncode}).\n"
        f"stderr: {result.stderr[:2000]}")
    assert "could not be located" not in output.lower(), (
        "exhaustive_find_cell failed.")
    assert "lost particle" not in output.lower(), (
        "Lost particle during transport.")

    # 阶段 2: 几何初始化
    assert "Entering cell" in output, (
        "No 'Entering cell' — find_cell_inner verbose missing.")

    # 阶段 3: 相变事件
    entering_events = output.count("Entering stochastic media")
    exiting_events = output.count("Exiting stochastic media")

    assert entering_events > 0, (
        f"No matrix→particle transitions detected ({entering_events}).")
    assert exiting_events > 0, (
        f"No particle→matrix transitions detected ({exiting_events}).")
    # 允许少量不对称：粒子在代际间被杀死/分裂，
    # 或模拟结束时仍在随机介质内部，因此不完全配对是正常的。
    # 2% 的容差足以排除结构性泄漏，同时兼容 PRNG 序列偏移。
    total = entering_events + exiting_events
    asymmetry = abs(entering_events - exiting_events) / total
    assert asymmetry < 0.02, (
        f"Phase asymmetry too large: {entering_events} entering vs "
        f"{exiting_events} exiting (ratio={asymmetry:.6f})")

    # 阶段 4: 表面穿越
    assert "Crossing surface" in output or "Reflected from surface" in output, (
        "No surface crossing events.")

    # 阶段 5: StatePoint
    sp_files = list(tmp_path.glob('statepoint.*.h5'))
    assert len(sp_files) == 1, (
        f"Expected 1 statepoint, found {len(sp_files)}")
    sp = openmc.StatePoint(sp_files[0])
    assert sp.keff.n > 0, f"keff should be positive, got {sp.keff}"
    assert 0.01 < sp.keff.n < 5.0, (
        f"keff {sp.keff.n} outside reasonable range")
