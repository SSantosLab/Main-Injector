from setuptools import setup, find_packages


setup(
    name='MainInjectr',
    version='4.0',
    description="MainInjector tools to Search for GW Counterparts",
    author="SSantosLab",
    packages=find_packages(),
    include_package_data=True,
    entry_points = {
        'main-injector scripts': [
            'main-injector run = main-injector.test:hello_world'
        ]
    }
)