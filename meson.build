project(
  'ad9361',
  'c',
  version: '0.3',
  default_options: [
    'c_std=c99',
    'warning_level=3',
  ],
)

add_project_arguments(
  [
    '-D_GNU_SOURCE=1',
    '-D_POSIX_C_SOURCE=200809L',
    '-DLIBAD9361_EXPORTS=1',
  ],
  language: 'c',
)

cc = meson.get_compiler('c')

libiio_dep = dependency('libiio')
m_dep = cc.find_library('m', required: false)

inc = include_directories('.')

ad9361_lib = library(
  'ad9361',
  [
    'ad9361.h',
    'ad9361_multichip_sync.c',
    'ad9361_baseband_auto_rate.c',
    'ad9361_design_taps.c',
    'ad9361_calculate_rf_clock_chain.c',
    'ad9361_fmcomms5_phase_sync.c',

    'filterdesigner/internal_design_filter_cg.c',
    'filterdesigner/internal_design_filter_cg_emxutil.c',
    'filterdesigner/rtGetInf.c',
    'filterdesigner/rtGetNaN.c',
    'filterdesigner/rt_nonfinite.c',
  ],
  include_directories: inc,
  dependencies: [libiio_dep, m_dep],
  install: true,
)

subdir('test')
