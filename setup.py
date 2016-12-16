from setuptools import setup

setup(
  name='moltemplate',
  packages=['moltemplate'],
  version='2.0.4',
  description='A general cross-platform text-based molecule builder for LAMMPS',
  author='Andrew Jewett',
  author_email='jewett.aij@gmail.com',
  url='https://github.com/jewettaij/moltemplate',
  download_url='https://github.com/jewettaij/moltemplate/tarball/v2.0.4',
  keywords=['simulation', 'lammps', 'molecule', 'builder'],
  license='BSD 3-Clause License',
  scripts=['moltemplate/scripts/moltemplate.sh'],
  entry_points={
    'console_scripts': [
        'ttree.py=moltemplate.ttree:main',
        'ttree_render.py=moltemplate.ttree_render:main',
        'lttree.py=moltemplate.lttree:main',
        'lttree_check.py=moltemplate.lttree_check:main',
        'lttree_postprocess.py=moltemplate.lttree_postprocess:main',
        'nbody_by_type.py=moltemplate.nbody_by_type:main',
        'nbody_fix_ttree_assignments.py=moltemplate.nbody_fix_ttree_assignments:main',
        'nbody_reorder_atoms.py=moltemplate.nbody_reorder_atoms:main',
        'pdbsort.py=moltemplate.pdbsort:main',
        'postprocess_input_script.py=moltemplate.postprocess_input_script:main',
        'remove_duplicate_atoms.py=moltemplate.remove_duplicate_atoms:main',
        'remove_duplicates_nbody.py=moltemplate.remove_duplicates_nbody:main',
        'renumber_DATA_first_column.py=moltemplate.renumber_DATA_first_column:main',
        'dump2data.py=moltemplate.dump2data:main',
        'raw2data.py=moltemplate.raw2data:main',
        'bonds_by_type.py=moltemplate.bonds_by_type:main']},
  package_data={'moltemplate': ['force_fields/*.lt']},
  # install_requires=['numpy', 'scipy', 'biopython'],
  setup_requires=['pytest-runner'],
  tests_require=['pytest', 'pandas'],
  zip_safe=True,
  include_package_data=True
)
