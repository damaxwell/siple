from distutils.core import setup


setup(
    name='siple',
    description='siple: a small inverse problems library',
    version='0.1',
    packages=['siple', 'siple.gradient', 'siple.linalg', 'siple.opt'],
    license='GNU Public License v2',
    author="David Maxwell",
    author_email="damaxwell@alaska.edu",
    long_description=open('README').read(),
)
