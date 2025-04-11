from setuptools import setup,find_packages
setup(name='bioencoder',
      version='1.0.0',
      description='A tool to encode AA to calculable encodings',
      author='iceshadows',
      author_email='adron_iceshadows@hotmail.com',
      zip_safe=False,
      packages=['bioencoder'],
      include_package_data=True,
      package_data={
        "": [".data/*.txt"]
      },
      license="AGPL"
      )