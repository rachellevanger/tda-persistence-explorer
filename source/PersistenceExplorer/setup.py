from setuptools import setup

setup(name='PersistenceExplorer',
      version='0.1',
      description='Jupyter-Embeddable Persistence Tool',
      url='',
      author='Rachel Levanger and Shaun Harker',
      author_email='rachel.levanger@rutgers.edu, shaun.harker@rutgers.edu',
      license='MIT',
      packages=['PersistenceExplorer'],
      install_requires=[
          'IPython'
      ],
      scripts=['bin/ImagePersistence'],
      zip_safe=False)
