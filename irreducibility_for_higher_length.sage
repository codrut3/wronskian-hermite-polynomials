# Copyright 2020 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
###############################################################################

# This program checks irreducibility for the remainder polynomial for
# partitions of length L <= 5.

from sage.symbolic.expression_conversions import polynomial

def hermite_poly(n, x):
  H = hermite(n, x)
  H_rescaled = H(x/sqrt(2))/(2^(n/2))
  H_rescaled = H_rescaled.expand()
  return polynomial(H_rescaled, base_ring=QQ)

def degree_vector(p):
  return [part + p.length() - i - 1 for i, part in enumerate(p)]

def compute_delta(p):
  degrees = degree_vector(p)
  S = Subsets(degrees, 2)
  D = 1
  for n_i, n_j in S:
    D *= abs(n_j - n_i)
  return D

def core_size(p):
  even = 0
  odd = 0
  degrees = degree_vector(p)
  for n_i in degrees:
    if n_i%2 == 0:
      even += 1
    else:
      odd += 1
  return ((odd - even) * (odd - even + 1))/2

x = PolynomialRing(QQ, 'x').gen()
R.<x> = ZZ[]

def compute_wronskian_hermite(p):
  if p.length() == 1:
    return hermite_poly(p.size(), x)
  else:
    degrees = degree_vector(p)
    degrees.reverse()

    W = wronskian(*[hermite_poly(n, x) for n in degrees])
    H_lambda = W/compute_delta(p)
    return H_lambda


# Pairs (L, N) for which to check irreducibility.
# All partitions of n <= N of length L will be tested.
to_check = [(2, 1000), (3, 250), (4, 150), (5, 110)]
for L, N in to_check:
  for n in list(range(1, N+1)):
    print('Checking partition size %s for length %s' % (n, L))
    partition_list = Partitions(n, length=L).list()
    for p in partition_list:
      R_lambda = compute_wronskian_hermite(p)
      s = core_size(p)
      R_lambda = R_lambda / x^s
      R_lambda = R(R_lambda)
      if not R_lambda.is_irreducible():
        print(p, s, R_lambda)
