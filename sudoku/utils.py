"""
Utility classes/functions for Sudoku solving.

Abhinav Sinha
2019-12-03
"""
from enum import Enum

from backtracking import SudokuBacktracker

class Algorithms(Enum):
    backtrack = 1
    constraint = 2


def read_puzzle(alg):
    """
    Read puzzle from STDIN. The expected format is: 
        row,col,value
    for each given in the puzzle. All unspecified cells are assumed to be blank.
    When all givens have been entered, input: done.
    """
    print("Order?", end = " ")
    n = int(input())
    sl = n**2 # side length
    print("Verify board has side length:", sl)
    print("Enter givens as\n\trow,col,value\nfollowed by\n\tdone\nwhen complete:")

    givens = dict()
    s = input()
    while s != "done":
        s = s.split(",")
        c = int(s[0]) * sl + int(s[1])
        givens[c] = int(s[2])
        s = input()

    if alg == Algorithms.backtrack:
        return SudokuBacktracker(n, givens)
    elif alg == Algorithms.constraint:
        # TODO
        pass

