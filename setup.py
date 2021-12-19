from setuptools import setup, find_packages

setup(
    name='glide',
    version='1',
    description='Denoising graphs using diffusion state distance link prediction',
    author='Henri Schmidt and Kapil Devkota',
    author_email='henri.schmidt@tufts.edu',
    url='https://github.com/kap-devkota/Trimming_Functional',
    packages=find_packages(exclude=('tests', 'docs', 'results', 'data')),
    package_dir={'glide':'glide'},
    package_data={'glide':['data/*.dat']}
)
