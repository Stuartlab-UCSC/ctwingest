from setuptools import setup, find_packages

setup(
    name='ctwingest',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "scanpy==1.3.7", "scipy", "pandas",
        "beautifulsoup4", "requests", "requests_toolbelt", "urllib3==1.26.5"
    ],
)
