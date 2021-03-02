"""
import click
import sys
from connect_the_dots.io import get_CZI_metadata


@click.command()
@click.option('--path')#, help='Path to the file of interest.')
@click.option('--name')#, help='Name of the the file')
def load_zstack(path, name):
    #Simple program that greets NAME for a total of COUNT times.
    #click.echo(f'{path} \n {name}!')
    print(f'{path} \n {name}!')
    
"""        
import click

@click.command()
@click.option('--path')
@click.option('--name')
def hello(path):
    print(f'Hello World! {path} {name}')
    #click.echo(f'Hello World! {path} {name}')
    #click.echo(path)
    #click.echo(name)