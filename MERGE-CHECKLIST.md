*Code Quality Checks*

- Abstraction Go down one abstraction level per function call
- Multi step processes should be constructed in a way that makes it impossible to skip a step. For example, a read operation should only be available as a method of an object representing a successfully  opened-file.
- Code has no side effects
- only one exit point from a function
- Initialize all variables
- Wrap calls to libraries for easier replacement
- All names in code (variables, functions, classes etc) are meaningful and used for what they are named for
- Function length is maximum 20 lines
- No more than one loop in a function
- Functions have maximum of three parameters
- Classes have few members and methods.
- All functions that can fail report the failure and the caller handles the error
- The code assumes that everything that can fail, will fail.
- All  input is checked for correctness
- comments are on what, not how.
- Anything the can be const , make const.
- All code is covered with unit tests, (and tests pass) - 100% of code lines.

- Error messages contain:
  - What the software tried to do
  - What went wrong
  - How to remedy the problem
  - As many relevant details as are at hand at the time the error occurred
