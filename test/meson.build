test(
  'FilterDesigner works',
  executable(
    'filter_designer_test',
    'filter_designer_test.c',
    link_with: [ad9361_lib],
    dependencies: [libiio_dep],
    include_directories: inc,
  ),
  workdir: meson.current_source_dir(),
)

test(
  'Generating rates works',
  executable(
    'generate_rates_test',
    'gen_rates_test.c',
    link_with: [ad9361_lib],
    dependencies: [libiio_dep],
    include_directories: inc,
  ),
  timeout: 60,
)

test(
  'FilterDesigner hardware works',
  executable(
    'filter_designer_hw_test',
    'filter_designer_hw.c',
    link_with: [ad9361_lib],
    dependencies: [libiio_dep],
    include_directories: inc,
  ),
)

test(
  'FMComms5 syncing works',
  executable(
    'fmcomms5_sync_test',
    'fmcomms5_sync_test.c',
    link_with: [ad9361_lib],
    dependencies: [libiio_dep, m_dep],
    include_directories: inc,
  ),
)

test(
  'Auto rate works',
  executable(
    'auto_rate_test',
    'auto_rate_test_hw.c',
    link_with: [ad9361_lib],
    dependencies: [libiio_dep],
    include_directories: inc,
  ),
)
