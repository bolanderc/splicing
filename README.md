# splicing

This package provides a tool to splice two near-field signatures together. This allows for a low-fidelity aerodynamic tool to be used in conjunction with a high-fidelity tool to make small changes to a geometry in a given region. Oftentimes there are regions of supersonic aircraft that are unsuitable for low-fidelity tools and this tool allows for portions of the aircraft geometry to be used for small changes.

```python
import splicing
import numpy as np

front_sig = np.genfromtxt('frontsig')
rear_sig = np.genfromtxt('rearsig')
x_cut = 0.5
blending = [10, 'linear']
l_ref = 32.92
splice = splicing.Splicing(front_sig, rear_sig, x_cut, blending, l_ref)
spliced_nearfield = splice.splice_sigs()

```

## Notes

splicing was supported by the NASA University Leadership Initiative (ULI) program under federal award number NNX17AJ96A, titled Adaptive Aerostructures for Revolutionary Civil Supersonic Transportation.

## Documentation

See doc strings in code.

## Installation

You can either download the source as a ZIP file and extract the contents or clone the pyldb repository using Git.

### Downloading source as a ZIP file

1. Open a web browser and navigate to <https://github.com/bolanderc/splicing
2. Make sure the branch is set to 'Master'
3. Click the `Clone or download` button
4. Select `Download ZIP`
5. Extract the downloaded ZIP file to a local directory on your machine

### Cloning the Github repository

1. From the command prompt, navigate to the directory where pyldb will be installed
2. `git clone https://github.com/bolanderc/splicing

## Testing

Unit tests are implemented using the pytest module and are run using the following command in the base directory.

```
pytest test.py
```

## Support

Contact [christian.bolander@aggiemail.usu.edu](mailto:christian.bolander@aggiemail.usu.edu) with any questions.

## License

This project is licensed under the MIT license. See LICENSE file for more information.