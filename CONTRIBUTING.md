# Contributing

These are guidelines for contributing to the RTMsim module. Contributors can open issues or create pull requests.

## Issues

Contributors can open new issues on the GitHub repository [here](https://github.com/obertscheiderfhwn/RTMsim/issues).
An issue can either report an existing bug, or request a new feature.

### Bug Reports

A new bug report can be opened [here](https://github.com/obertscheiderfhwn/RTMsim/issues/new?template=bug_report.md).

### Feature Requests

A new feature request can be opened [here](https://github.com/obertscheiderfhwn/RTMsim/issues/new?template=feature_request.md).
Feature requests that are deemed feasible will be considered by the developers, and could even be addressed by a contributor through a pull request.
Feature requests that are deemed infeasible will likely be denied.

## Pull Requests

Contributors can propose changes to the code in the repository by creating a pull request as follows:

- Fork the base repository [here](https://github.com/obertscheiderfhwn/RTMsim/fork).
- Clone the forked repository, make changes, and push them back to the fork.
- Create a pull request between the base and forked repositories [here](https://github.com/obertscheiderfhwn/RTMsim/pulls).
- Wait for the pull request to be either approved or dismissed. Approval and subsequent merging of pull requests is contingent upon:

  - Existing tests are successful and new tests are added and successful.
  - The changes are properly documented.
  - The changes provide an appropriate and substantial improvement to the repository.

## Testing
Please verify that the code still produced the expected results for the basic V&V cases. The validation and verification cases can be executed in the GUI with the input files input_case1_coarsemesh.txt, input_case1_finemesh.txt, input_case2_coarsemesh.txt, input_case2_finemesh.txt, input_case3_coarsemesh.txt, input_case3_finemesh.txt. The input file names are saved in directory inputfiles.
