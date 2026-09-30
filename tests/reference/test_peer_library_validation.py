"""Executable validation of independent KerrGeoPy and PyGRO oracles."""

from __future__ import annotations

from importlib import metadata

import numpy as np
import pytest

from .peer_oracles import (
    case_by_id,
    exact_schwarzschild_circular_reference,
    integrate_with_pygro,
    kerrgeopy_bound_reference,
    kerrgeopy_circular_reference,
    load_cases,
    make_pygro_engine,
)


@pytest.fixture(scope="module")
def pygro_backend():
    pytest.importorskip("pygro", reason="optional peer library PyGRO is not installed")
    return make_pygro_engine()


def test_fixture_records_explicit_conventions_and_pending_review():
    manifest = load_cases()
    assert manifest["review_status"] == "pending"
    assert manifest["conventions"]["parameter"] == "Timelike proper time tau"
    assert {item["record_id"] for item in manifest["provenance"]} == {
        "INF-006",
        "INF-007",
    }
    assert all(item["manual_review"] == "pending" for item in manifest["provenance"])


def test_installed_peer_versions_match_manifest():
    pytest.importorskip("kerrgeopy", reason="optional peer library KerrGeoPy is not installed")
    pytest.importorskip("pygro", reason="optional peer library PyGRO is not installed")
    manifest = load_cases()["libraries"]
    assert metadata.version("kerrgeopy") == manifest["kerrgeopy"]["tested_version"]
    assert tuple(map(int, metadata.version("pygro").split("."))) >= tuple(
        map(int, manifest["pygro"]["minimum_version"].split("."))
    )


def test_stable_schwarzschild_circular_shared_state(pygro_backend):
    pytest.importorskip("kerrgeopy", reason="optional peer library KerrGeoPy is not installed")
    case = case_by_id("schwarzschild_circular_r10")
    exact = exact_schwarzschild_circular_reference(case["radius_over_M"])
    reference = kerrgeopy_circular_reference(case)
    np.testing.assert_allclose(reference.coordinates, exact.coordinates, rtol=0.0, atol=2e-12)
    np.testing.assert_allclose(reference.four_velocity, exact.four_velocity, rtol=0.0, atol=2e-12)

    metric, engine = pygro_backend
    result = integrate_with_pygro(
        metric,
        engine,
        reference,
        spin=case["a_over_M"],
        goal=case["pygro_goal"],
    )
    assert result.max_abs_endpoint_error <= case["tolerances"]["coordinate_max_abs"]
    assert (
        result.max_abs_trajectory_error
        <= case["tolerances"]["trajectory_coordinate_max_abs"]
    )
    assert result.max_abs_normalization_error <= case["tolerances"]["normalization_max_abs"]


def test_schwarzschild_isco_against_exact_circular_state(pygro_backend):
    case = case_by_id("schwarzschild_isco_r6")
    assert case["kerrgeopy_applicable"] is False
    reference = exact_schwarzschild_circular_reference(case["radius_over_M"])
    metric, engine = pygro_backend
    result = integrate_with_pygro(
        metric,
        engine,
        reference,
        spin=case["a_over_M"],
        goal=case["pygro_goal"],
    )
    radius_error = float(np.max(np.abs(result.coordinates[:, 1] - case["radius_over_M"])))
    assert result.max_abs_endpoint_error <= case["tolerances"]["coordinate_max_abs"]
    assert (
        result.max_abs_trajectory_error
        <= case["tolerances"]["trajectory_coordinate_max_abs"]
    )
    assert radius_error <= case["tolerances"]["radius_max_abs"]
    assert result.max_abs_normalization_error <= case["tolerances"]["normalization_max_abs"]


def test_eccentric_schwarzschild_one_radial_period_and_precession(pygro_backend):
    pytest.importorskip("kerrgeopy", reason="optional peer library KerrGeoPy is not installed")
    pytest.importorskip("scipy", reason="SciPy is required for Mino-to-proper-time mapping")
    case = case_by_id("schwarzschild_eccentric")
    reference = kerrgeopy_bound_reference(case)
    metric, engine = pygro_backend
    result = integrate_with_pygro(
        metric,
        engine,
        reference,
        spin=case["a_over_M"],
        goal=case["pygro_goal"],
    )
    delta_phi = float(reference.coordinates[3, -1] - reference.coordinates[3, 0])
    periapsis_advance = delta_phi - 2.0 * np.pi
    assert periapsis_advance >= case["tolerances"]["minimum_periapsis_advance_rad"]
    assert result.max_abs_endpoint_error <= case["tolerances"]["coordinate_max_abs"]
    assert (
        result.max_abs_trajectory_error
        <= case["tolerances"]["trajectory_coordinate_max_abs"]
    )
    assert result.max_abs_normalization_error <= case["tolerances"]["normalization_max_abs"]


def test_bound_kerr_shared_state_and_tolerance_convergence(pygro_backend):
    pytest.importorskip("kerrgeopy", reason="optional peer library KerrGeoPy is not installed")
    pytest.importorskip("scipy", reason="SciPy is required for Mino-to-proper-time mapping")
    case = case_by_id("kerr_stable_bound")
    reference = kerrgeopy_bound_reference(case)
    metric, engine = pygro_backend
    coarse_goal, fine_goal = case["convergence_goals"]
    coarse = integrate_with_pygro(
        metric,
        engine,
        reference,
        spin=case["a_over_M"],
        goal=coarse_goal,
    )
    fine = integrate_with_pygro(
        metric,
        engine,
        reference,
        spin=case["a_over_M"],
        goal=fine_goal,
    )
    assert fine.max_abs_endpoint_error <= case["tolerances"]["coordinate_max_abs"]
    assert (
        fine.max_abs_trajectory_error
        <= case["tolerances"]["trajectory_coordinate_max_abs"]
    )
    assert fine.max_abs_normalization_error <= case["tolerances"]["normalization_max_abs"]
    assert fine.max_abs_endpoint_error < coarse.max_abs_endpoint_error
    assert (
        fine.max_abs_endpoint_error / coarse.max_abs_endpoint_error
        <= case["tolerances"]["fine_to_coarse_error_ratio_max"]
    )
    assert fine.max_abs_normalization_error < coarse.max_abs_normalization_error
