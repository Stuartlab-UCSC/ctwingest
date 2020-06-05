from setuptools import setup, find_packages

setup(
    name='ctwingest',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "anndata==0.7.1", "scanpy==1.4.4", "scipy", "pandas",
        "beautifulsoup4", "requests", "requests_toolbelt", "urllib3==1.24.2"
    ],
)
