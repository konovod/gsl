# GNU GSL — fork

[![CMake](https://github.com/konovod/gsl/actions/workflows/cmake.yml/badge.svg)](https://github.com/konovod/gsl/actions/workflows/cmake.yml)
[![Docs](https://github.com/konovod/gsl/actions/workflows/docs.yml/badge.svg)](https://konovod.github.io/gsl/)

Contrary to everything else in this fork, this file is written by hand.
This is a fork of GSL.
I dreamed for a long time to make a GSL well-maintained and modern library. And now i'm kind of fullfilling it with a vibecoding.

## Contra:

Believe me, i know how horrible it sounds and yes it IS full of AI-slop with minimal review of a random guy who make important design decisions by clicking to Recommended solution.
It is even done not with a top-tier model (most of code written with deepseek-4.1-flash), so expect not just unmaintainable code but also a factual bugs.

## Pro (what is done in this fork):

- CMake support to build it in a modern way on all supported platform (Windows/Linux/macOS)
- Github Actions to actually test in on all supported platforms
- Right now i'm going to fix all 219 open issues of GSL bug tracker. All reproducible bugs are already fixed, now in process of checking performance related and new feature requests. I'm trying to make fixes a separate commits to be able at least in theory merge them to upstream.

## What to look

Agent keeps list of all changes in `FORKNEWS`.
Also there is `SAVANNA_REVIEW.md` and `SAVANNA_TRIAGE.md` - this is AI-generated summaries of work on issues in tracker.
There is an AI-generated `CMake.md` about building process, you can also check github actions to understand how to build it on your system.
