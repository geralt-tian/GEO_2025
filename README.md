# Efficient and High-Accuracy Secure Two-Party Protocols for a Class of Functions with Real-number Inputs



## Setup

For setup instructions, please refer to the README file located in the `SCI` folder.

We successfully completed the compilation on Ubuntu 22.04.5 LTS with Intel(R) Xeon(R) Platinum.


## Code Structure

The project is organized as follows:

- **/SCI/tests**
  Contains all our related code, including implementations of activation functions and models.

- **/SCI/tests/activation**
  The **/SCI/tests/activation** directory includes our evaluation functions:
  - MW
  - exp
  - sin
  - division
  - softmax

## Running Tests

To run the unit tests, use the following command:

```bash
./SCI/build/bin/our-MW r=1 & ./SCI/build/bin/our-MW r=2
```

**Reference Repository:**  
**Project webpage:** <[BOLT](https://github.com/Clive2312/EzPC/tree/bert/SCI)>

**Reference Repository:**  
**Project webpage:** <[SEAF](https://github.com/geralt-tian/SEAF)>