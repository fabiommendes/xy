import sys

from invoke import task

python = sys.executable


@task
def build(ctx):
    ctx.run(f'{python} setup.py build_ext --inplace')


@task
def test(ctx):
    build(ctx)
    ctx.run(f'{python} -m pytest')


@task
def coverage(ctx):
    build(ctx)
    ctx.run(f'{python} -m pytest --cov')
