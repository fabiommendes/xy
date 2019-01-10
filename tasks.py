import os
import sys

from invoke import task

python = sys.executable


@task
def build(ctx, disable_cython=False):
    env = dict(os.environ)
    if disable_cython:
        env['DISABLE_CYTHON'] = 'true'
    ctx.run(f'{python} setup.py build_ext --inplace', env=env)


@task
def test(ctx):
    build(ctx)
    ctx.run(f'{python} -m pytest')


@task
def coverage(ctx):
    build(ctx)
    ctx.run(f'{python} -m pytest --cov')
