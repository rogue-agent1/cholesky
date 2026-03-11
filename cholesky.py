#!/usr/bin/env python3
"""Cholesky decomposition — A = LLᵀ for positive definite matrices."""
import sys, math
def cholesky(A):
    n = len(A); L = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            s = sum(L[i][k]*L[j][k] for k in range(j))
            L[i][j] = math.sqrt(A[i][i] - s) if i == j else (A[i][j] - s) / L[j][j]
    return L
if __name__ == "__main__":
    A = [[4,12,-16],[12,37,-43],[-16,-43,98]]
    L = cholesky(A)
    print("A:"); [print(f"  {r}") for r in A]
    print("L:"); [print(f"  {[round(x,3) for x in r]}") for r in L]
    # Verify L·Lᵀ = A
    n = len(L)
    LLt = [[sum(L[i][k]*L[j][k] for k in range(n)) for j in range(n)] for i in range(n)]
    print("LLᵀ:"); [print(f"  {[round(x,3) for x in r]}") for r in LLt]
