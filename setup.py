from setuptools import setup, find_packages


config = {
    'name': 'MainInjector',
    'version': '4.0',
    'description': 'MainInjector tools to Search for GW Counterparts',
    'author': 'SSantosLab',
    'packages': find_packages(),
    'include_package_data': True,
    'classifiers': [
        "Programming Language :: Python :: 3",
        "Programming Language :: Fortran :: 77",
        "Programming Language :: C",
        'License :: ',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics'
    ],
    'entry_points': {
        'console_scripts' : ['main-injector-config=main_injector.cli.make_config:main']
    }
}
setup(**config)