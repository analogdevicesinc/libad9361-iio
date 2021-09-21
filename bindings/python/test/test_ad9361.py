import pytest

import iio
import ad9361
import os
import pathlib

ad9361_all = ["fmcomms2", "fmcomms5", "ad9361", "ad9364"]
path = pathlib.Path(__file__).parent.absolute()
filterfile_name = os.path.join(path, "LTE15_MHz.ftr")


@pytest.mark.iio_hardware(ad9361_all)
@pytest.mark.parametrize("rate_msps", [*(range(1, 30))])
def test_bb_rate(iio_uri, rate_msps):
    rate = rate_msps * 1e6
    ctx = iio.Context(iio_uri)
    dev = ctx.find_device("ad9361-phy")
    ad9361.set_bb_rate(dev, int(rate))
    hw = float(dev.find_channel("voltage0").attrs["sampling_frequency"].value)
    assert abs(rate - hw) < 3


@pytest.mark.iio_hardware(ad9361_all)
@pytest.mark.parametrize("rate_msps", [*(range(1, 30))])
def test_bb_rate_custom_filter_manual(iio_uri, rate_msps):
    rate = rate_msps * 1e6
    ctx = iio.Context(iio_uri)
    dev = ctx.find_device("ad9361-phy")
    Fstop = rate * 0.6
    Fpass = rate * 0.4
    wnom = Fpass
    ad9361.set_bb_rate_custom_filter_manual(
        dev, int(rate), int(Fpass), int(Fstop), int(wnom), int(wnom)
    )
    hw = float(dev.find_channel("voltage0").attrs["sampling_frequency"].value)
    assert abs(rate - hw) < 3


@pytest.mark.iio_hardware(ad9361_all)
@pytest.mark.parametrize("rate_msps", [*(range(1, 30))])
def test_bb_rate_custom_filter(iio_uri, rate_msps):
    rate = rate_msps * 1e6
    ctx = iio.Context(iio_uri)
    dev = ctx.find_device("ad9361-phy")
    ad9361.set_bb_rate_custom_filter_auto(dev, int(rate))
    hw = float(dev.find_channel("voltage0").attrs["sampling_frequency"].value)
    assert abs(rate - hw) < 3


@pytest.mark.iio_hardware("fmcomms5")
def test_mcs(iio_uri):
    ctx = iio.Context(iio_uri)
    main = ctx.find_device("ad9361-phy")
    assert main
    secondary = ctx.find_device("ad9361-phy-B")
    assert secondary
    ret = ad9361.multichip_sync(main, [secondary], 3)
    assert ret == 0


@pytest.mark.iio_hardware("fmcomms5")
def test_fmc5_mcs(iio_uri):
    ctx = iio.Context(iio_uri)
    ret = ad9361.fmcomms5_multichip_sync(ctx, 3)
    assert ret == 0


@pytest.mark.iio_hardware(ad9361_all)
def test_trx_fir_enable(iio_uri):
    ctx = iio.Context(iio_uri)
    dev = ctx.find_device("ad9361-phy")
    assert dev

    # Set sample rate so filter isn't needed
    dev.find_channel("voltage0").attrs["sampling_frequency"].value = "3000000"

    # disable filter and check
    ad9361.set_trx_fir_enable(dev, 0)
    assert dev.find_channel("voltage0").attrs["filter_fir_en"].value == "0"
    assert dev.find_channel("voltage0", True).attrs["filter_fir_en"].value == "0"

    # Load filter
    with open(filterfile_name, "r") as file:
        data = file.read()
    dev.attrs["filter_fir_config"].value = data

    # Enable filter and check
    ad9361.set_trx_fir_enable(dev, 1)
    assert ad9361.get_trx_fir_enable(dev) == 1
    assert dev.find_channel("voltage0").attrs["filter_fir_en"].value == "1"
    assert dev.find_channel("voltage0", True).attrs["filter_fir_en"].value == "1"

    # Disable filter and check
    ad9361.set_trx_fir_enable(dev, 0)
    assert dev.find_channel("voltage0").attrs["filter_fir_en"].value == "0"
    assert dev.find_channel("voltage0", True).attrs["filter_fir_en"].value == "0"


@pytest.mark.iio_hardware("fmcomms5")
def test_fmc5_phase_sync(iio_uri):
    ctx = iio.Context(iio_uri)
    ret = ad9361.fmcomms5_phase_sync(ctx, 1e9)
    assert ret == 0


def test_generate_fir_taps():

    rate = 10e6

    [rxfdp, txfdp] = ad9361.calculate_rf_clock_chain_fdp(rate)

    assert rxfdp.Rdata == rate
    assert txfdp.Rdata == rate

    [tapsRX, gain] = ad9361.generate_fir_taps(rxfdp)
    [tapsTX, gainTX] = ad9361.generate_fir_taps(txfdp)

    def compare(taps):
        rtaps = taps[int(len(taps) / 2) :]
        rtaps = rtaps[::-1]
        for k in range(64):
            assert taps[k] == rtaps[k]

    compare(tapsRX)
    compare(tapsTX)


def test_generate_clock_chain():

    rate = 10e6
    rate_gov = 0

    [rx_path_clks, tx_path_clks] = ad9361.calculate_rf_clock_chain(rate, rate_gov)

    rx_ref = [960000000, 480000000, 160000000, 80000000, 40000000, 10000000]
    tx_ref = [960000000, 240000000, 80000000, 40000000, 40000000, 10000000]
    assert rx_path_clks == rx_ref
    assert tx_path_clks == tx_ref
