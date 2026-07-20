==============
Release policy
==============

surface-pd is currently a pre-1.0 Beta package. It is used for active research
workflows, while its public interface and documentation continue to mature.

Pre-1.0 releases may introduce breaking API changes. Such changes must be
identified in release notes, accompanied by migration guidance, and reflected
in the package version before a release is published. Once the project is ready
to promise a stable public API, it will adopt a 1.0 or later version and update
its development-status metadata accordingly.

All user-visible version reporting is derived from the authoritative version in
``surface_pd/_version.py``. Project metadata, Python runtime metadata, command
line output, built wheels, and the Sphinx documentation must agree with it.
