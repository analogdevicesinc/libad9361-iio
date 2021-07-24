import pytest

import iio
import ad9361


@pytest.mark.parametrize("rate_msps", [*(range(1, 30))])
def test_bb_rate(rate_msps):
    rate = rate_msps * 1e6
    ctx = iio.Context("ip:192.168.86.39")
    dev = ctx.find_device("ad9361-phy")
    ad9361.set_bb_rate(dev, int(rate))
    hw = float(dev.find_channel("voltage0").attrs["sampling_frequency"].value)
    assert abs(rate - hw) < 3
