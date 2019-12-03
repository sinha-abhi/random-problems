"""
Sudoku solving via backtracing (aka brute force).

Abhinav Sinha
2019-12-03
"""
import time

import utils

from typing import Mapping


class SudokuBacktracer:
    def __init__(self, order, givens: Mapping[int, int]):
        """
        order   order of puzzle
                A puzzle of order n forms a n^2 by n^2 grid with n^4 values,
                each value being between 1 to n^2. The most common Sudoku is
                of order 3.
        givens  mapping of cell number to its value
                The cell numbers of a sudoku of order 3 are as follows:
                1   2   ...   9
                10  11  ...   18
                 .      .     .
                 .       .    .
                 .        .   .
                73  74  ...   81

        """
        self.givens = givens
        self.iter = 0
        self.order = order
        self.slen = order**2
        # init empty board
        self.board = [[0] * self.slen for i in range(self.slen)]

    def __backtrack(self):
        self.iter += 1
        row, col = self.__find_empty()
        if (row == -1) and (col == -1): # board is non-empty
            return True

        for i in range(self.slen):
            _val = i + 1 # tentative value
            if self.__valid(row, col, _val):
                self.board[row][col] = _val
                if self.__backtrack():
                    return True

                self.board[row][col] = 0

        return False

    def __cage(self, row, col):
        """
        Finds the cage containing the given cell.
        """
        r, c = self.__cage_location(row, col)
        cage = list()
        for i in range(r, r + self.order):
            cage.extend(self.board[i][c : c + self.order])

        return cage

    def __column(self, col):
        return [row[col] for row in self.board]

    def __cage_location(self, row, col):
        return self.__cage_start(row), self.__cage_start(col)

    def __cage_start(self, num):
        return num // self.order * self.order

    def __find_empty(self):
        """
        Finds then next cell that has no assigned value.
        """
        for row in range(self.slen):
            for col in range(self.slen):
                if self.board[row][col] == 0:
                    return row, col

        return -1, -1

    def __init_board(self, verbose=False):
        for k, v in self.givens.items():
            row = k // self.slen
            col = k % self.slen
            self.board[row][col] = v
        if verbose:
            self.print_board()

    def __row(self, row):
        return self.board[row]

    def __valid(self, r, c, val):
        """
        Checks if the value is valid in the proposed position.
        """
        col = self.__column(c)
        row = self.__row(r)
        cage = self.__cage(r, c)

        return ((val not in row) and (val not in col) and (val not in cage))

    def print_board(self):
        print(*self.board, sep = "\n")

    def solve(self, verbose=False):
        print("Solving...")
        self.__init_board(verbose)
        pfs = time.perf_counter()
        pcs = time.process_time()
        self.__backtrack()
        pfe = time.perf_counter()
        pce = time.process_time()
        print("---------------------------")
        self.print_board()
        if verbose:
            print("---------------------------")
            print("Elapsed time:", pfe - pfs)
            print("CPU Process time:", pce - pcs)
            print("Iterations:", self.iter)
            print("---------------------------")


if __name__ == "__main__":
    puzzle = utils.read_puzzle(utils.Algorithms.backtrace)
    puzzle.solve(True)

