from __future__ import annotations

from pathlib import Path
import sys


HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from run_micro_batch_ablation import (  # noqa: E402
    build_schedules,
    fixed_batch_schedule,
    historical_jump_schedule,
)


def test_historical_schedule_matches_5_to_10_to_20_counts():
    schedule = historical_jump_schedule(
        seed_count=81,
        ten_percent_count=161,
        target_count=322,
    )
    assert schedule == [80, 161]
    assert 81 + sum(schedule) == 322


def test_fixed_micro_batches_stop_exactly_at_target():
    for batch_size in (32, 16, 8, 4):
        schedule = fixed_batch_schedule(
            seed_count=81,
            target_count=322,
            batch_size=batch_size,
        )
        assert 81 + sum(schedule) == 322
        assert max(schedule) <= batch_size
        assert all(value > 0 for value in schedule)


def test_preregistered_schedules_share_same_total_budget():
    schedules = build_schedules(
        seed_count=81,
        ten_percent_count=161,
        target_count=322,
    )
    assert set(schedules) == {
        "historical_jumps",
        "micro_32",
        "micro_16",
        "micro_8",
        "micro_4",
    }
    for schedule in schedules.values():
        assert 81 + sum(schedule) == 322


def test_smaller_micro_batches_require_more_refits():
    schedules = build_schedules(
        seed_count=81,
        ten_percent_count=161,
        target_count=322,
    )
    assert len(schedules["historical_jumps"]) == 2
    assert len(schedules["micro_32"]) < len(schedules["micro_16"])
    assert len(schedules["micro_16"]) < len(schedules["micro_8"])
    assert len(schedules["micro_8"]) < len(schedules["micro_4"])
