from setuptools import setup

setup(name='PersistenceExplorer',
      version='0.1',
      description='Jupyter-Embeddable Persistence Tool',
      url='',
      author='Rachel Levanger and Shaun Harker',
      author_email='flyingcircus@example.com',
      license='MIT',
      packages=['PersistenceExplorer'],
      install_requires=[
          'IPython'
      ],
      scripts=['bin/ImagePersistence'],
      zip_safe=False)
