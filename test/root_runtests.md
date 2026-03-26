# `test/runtests.jl`

This is the root test entrypoint for the repository-level Julia package.

## What it does

`test/runtests.jl`:

- loads `Test`
- includes the more detailed test definitions from the nested test file
- keeps the root package test command simple

## Why it matters

The root test file is the file `Pkg.test()` uses from the repository root. Keeping it small is intentional:

- it avoids duplicating test logic
- it keeps the package test path easy to follow
- it lets the actual test cases live in one place

## Relationship to the nested test file

The real assertions remain in the nested `BioToolkit.jl/test/runtests.jl` file. This root file just acts as the bridge so the repository-level package can run the same tests.

## Why this file is small

It does not need extra logic. Its only job is to connect the root package test command to the existing test suite.
