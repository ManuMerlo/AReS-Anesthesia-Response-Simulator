from setuptools import setup, find_packages
import os

# Ensure README.md exists before reading
readme_path = os.path.join(os.path.dirname(__file__), "README.md")
if not os.path.exists(readme_path):
    raise FileNotFoundError(f"README.md not found at {readme_path}")

# Load the README file for a detailed package description
with open(readme_path, "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='AReS',  # Package name
    version='0.1',  # Package version
    author='Manuela Merlo, Sara Hosseinirad',
    author_email='merlomanu1999@gmail.com, sarahrad@ece.ubc.ca',
    description='Open-source AReS for modeling patient responses under anesthesia.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/ManuMerlo/AReS-Patient-Simulator',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "pandas",
        "numpy",
        "control",
        "scipy"
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
    ],
    python_requires='>=3.9',
    entry_points={
        'console_scripts': [
            'ares-simulator=AReS.main:main',  # Allows running `ares-simulator` from terminal
        ],
    },
)
