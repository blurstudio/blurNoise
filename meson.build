project('blurNoise', 'cpp', default_options: ['cpp_std=c++20'])

fs = import('fs')
maya_dep = dependency('maya')
maya_name_suffix = maya_dep.get_variable('name_suffix')
maya_version = maya_dep.get_variable('maya_version')

if fs.is_file('src/version.h')
  message('Using existing version.h')
  version_h = 'src/version.h'
else
  git = find_program('git', native: true, required: true)
  version_h = vcs_tag(
    command: [git, 'describe', '--tags', '--match', 'v[0-9]*', '--dirty=+'],
    fallback: 'v0.0.1',
    input: 'src/version.h.in',
    output: 'version.h',
  )
endif

source_files = files([
  'src/OpenSimplex2F.cpp',
  'src/blurNoise.cpp',
  'src/pluginMain.cpp',
])

outlib = shared_library(
  meson.project_name(),
  source_files + version_h,
  install: true,
  install_dir : meson.global_source_root() / 'output_Maya' + maya_version,
  include_directories : include_directories(['src']),
  dependencies : maya_dep,
  name_prefix : '',
  name_suffix : maya_name_suffix,
)

