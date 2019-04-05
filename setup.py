from setuptools import setup

setup(

  name='moltemplate',

  packages=['moltemplate',
            'moltemplate.nbody_alt_symmetry'],

  package_dir={'moltemplate': 'moltemplate'},           #.py files are in "moltemplate/"

  package_data={'moltemplate': ['force_fields/*.lt']},  #.lt files are in "moltemplate/force_fields/"

  #package_data={'moltemplate/force_fields':['*.lt']}
  #
  #package_data={'moltemplate/force_fields':
  #              ['compass_published.lt',
  #               'cooke_deserno_lipid.lt',
  #               'gaff2.lt',
  #               'gaff.lt',
  #               'graphene.lt',
  #               'graphite.lt',
  #               'loplsaa.lt',
  #               'martini.lt',
  #               'oplsaa.lt',
  #               'sdk.lt',
  #               'spce_ice_rect16.lt',
  #               'spce_ice_rect32.lt',
  #               'spce_ice_rect8.lt',
  #               'spce.lt',
  #               'tip3p_1983_charmm.lt',
  #               'tip3p_1983.lt',
  #               'tip3p_2004.lt',
  #               'tip5p.lt',
  #               'trappe1998.lt',
  #               'watmw.lt']},

  description='A general cross-platform text-based molecule builder for LAMMPS',

  long_description='Moltemplate is a general cross-platform text-based molecule builder for LAMMPS and ESPResSo. Moltemplate was intended for building custom coarse-grained molecular models, but it can be used to prepare realistic all-atom simulations as well.  It supports a variety of force fields for all-atom and coarse-grained modeling (including many-body forces and non-point-like particles).  New force fields and examples are added continually by users.',

  author='Andrew Jewett',

  author_email='jewett.aij@gmail.com',

  url='https://github.com/jewettaij/moltemplate',

  download_url='https://github.com/jewettaij/moltemplate/archive/v2.10.12.zip',

  version='2.10.12',

  keywords=['simulation', 'LAMMPS', 'molecule editor', 'molecule builder',
            'ESPResSo'],
            

  # BSD 3-Clause License:
  # - http://choosealicense.com/licenses/bsd-3-clause
  # - http://opensource.org/licenses/BSD-3-Clause

  #license='BSD',

  classifiers=['Environment :: Console',
               'License :: OSI Approved :: BSD License',
               'Operating System :: MacOS :: MacOS X',
               'Operating System :: POSIX :: Linux',
               'Operating System :: Microsoft :: Windows',
               'Programming Language :: Python',
               'Programming Language :: Unix Shell',
               'Topic :: Scientific/Engineering :: Chemistry',
               'Topic :: Scientific/Engineering :: Physics',
               'Topic :: Multimedia :: Graphics :: 3D Modeling',
               'Intended Audience :: Science/Research'],

  scripts=['moltemplate/scripts/moltemplate.sh',
           'moltemplate/scripts/cleanup_moltemplate.sh',
           'moltemplate/scripts/emoltemplate.sh'],

  entry_points={
    'console_scripts': [
        'ttree.py=moltemplate.ttree:main',
        'ttree_render.py=moltemplate.ttree_render:main',
        'bonds_by_type.py=moltemplate.bonds_by_type:main',
        'charge_by_bond.py=moltemplate.charge_by_bond:main',
        'dump2data.py=moltemplate.dump2data:main',
        'extract_espresso_atom_types.py=moltemplate.extract_espresso_atom_types:main',
        'extract_lammps_data.py=moltemplate.extract_lammps_data:main',
        'ettree.py=moltemplate.ettree:main',
        'genpoly.py=moltemplate.ettree:main',
        'ltemplify.py=moltemplate.ltemplify:main',
        'lttree.py=moltemplate.lttree:main',
        'lttree_check.py=moltemplate.lttree_check:main',
        'lttree_postprocess.py=moltemplate.lttree_postprocess:main',
        'nbody_by_type.py=moltemplate.nbody_by_type:main',
        'nbody_fix_ttree_assignments.py=moltemplate.nbody_fix_ttree_assignments:main',
        'nbody_reorder_atoms.py=moltemplate.nbody_reorder_atoms:main',
        'pdbsort.py=moltemplate.pdbsort:main',
        'postprocess_input_script.py=moltemplate.postprocess_input_script:main',
        'postprocess_coeffs.py=moltemplate.postprocess_coeffs:main',
        'raw2data.py=moltemplate.raw2data:main',
        'remove_duplicate_atoms.py=moltemplate.remove_duplicate_atoms:main',
        'remove_duplicates_nbody.py=moltemplate.remove_duplicates_nbody:main',
        'renumber_DATA_first_column.py=moltemplate.renumber_DATA_first_column:main']},

  # install_requires=['numpy', 'scipy'],
  setup_requires=['pytest-runner'],
  tests_require=['pytest'],
  zip_safe=True,
  include_package_data=True
)
