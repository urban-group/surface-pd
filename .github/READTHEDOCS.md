# Read the Docs maintainer guide

The public documentation is hosted at
<https://surface-pd.readthedocs.io/en/latest/>. Read the Docs builds and
publishes the site from `urban-group/surface-pd`; GitHub Actions independently
checks the Sphinx HTML and doctest builders before publication.

## Standard configuration

- Keep the project connected to `https://github.com/urban-group/surface-pd`.
- Keep `latest` active and mapped to the `main` branch.
- Keep the repository-root `.readthedocs.yaml` as the build configuration.
- Use the standard Read the Docs GitHub integration. Do not create a duplicate
  RTD project, an API token deployment workflow, or a GitHub Actions deployment
  trigger.

## Verify automatic publication

After a push to `main`:

1. Open the RTD project and select **Builds**.
2. Confirm that exactly one `latest` build starts without a manual action.
3. Confirm that the build identifies the pushed commit from
   `urban-group/surface-pd` and completes successfully.
4. Open the public `latest` documentation and verify that it contains the
   pushed change.

GitHub Actions passing is necessary but does not prove publication. The RTD
build must also start and succeed for the same commit.

## Trigger a manual recovery build

Use a manual build to diagnose or republish the current `main` commit:

1. Open **Builds** and select **Build version**.
2. Choose `latest` and start the build.
3. Check that the clone step uses `urban-group/surface-pd` and that checkout
   resolves to the expected `main` commit.

A successful manual build proves that the repository and build configuration
work. It does not prove that push notifications work.

## Recover a missing push trigger

1. Open **Admin → Integrations** in the RTD project.
2. Open the GitHub integration and inspect **Recent Activity** for the push.
3. If there is no matching activity, re-sync the integration.
4. If re-syncing fails or the integration refers to another fork, remove that
   integration and create a new GitHub integration for
   `urban-group/surface-pd`.
5. Verify that the maintainer's RTD account is connected to GitHub and that the
   Read the Docs GitHub App has access to the organization repository. An
   organization owner may need to approve that access.
6. Trigger one manual `latest` build, then verify the integration with a later
   push to `main`.

If **Recent Activity** contains the push but no build starts, confirm that
`latest` is active and mapped to `main`. If a build starts but fails, diagnose
the first failing command in its build log before changing the integration.
