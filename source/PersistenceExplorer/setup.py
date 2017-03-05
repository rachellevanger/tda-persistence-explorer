
from setuptools import setup


setup(name='PersistenceExplorer',
      version='0.1',
      description='Jupyter-Embeddable Persistence Tool',
      url='https://github.com/rachellevanger/tda-persistence-explorer',
      author='Rachel Levanger and Shaun Harker',
      author_email='rachel@math.rutgers.edu, shaun.harker@rutgers.edu',
      license='MIT',
      packages=['PersistenceExplorer'],
      install_requires=[
          'IPython'
      ],
      include_package_data=True,
      data_files=[('',['bin/ImagePersistence'])], # includes binary executable
      zip_safe=False)
