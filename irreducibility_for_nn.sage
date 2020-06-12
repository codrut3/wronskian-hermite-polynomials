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

# This program checks irreducibility for the remainder polynomial for the
# partitions (n, n) where n <= 1000 or n is one of {2401, 4375}.

from sage.symbolic.expression_conversions import polynomial

def hermite_poly(n, x):
  H = hermite(n, x)
  H_rescaled = H(x/sqrt(2))/(2^(n/2))
  H_rescaled = H_rescaled.expand()
  return polynomial(H_rescaled, base_ring=QQ)

def degree_vector(p):
  return [part + p.length() - i - 1 for i, part in enumerate(p)]

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
    lambda_1 = p.get_part(0)
    lambda_2 = p.get_part(1)
    M = matrix(SR, 2, 2, [hermite_poly(lambda_2, x),
                          hermite_poly(lambda_1+1, x),
                          lambda_2*hermite_poly(lambda_2-1, x),
                          (lambda_1+1)*hermite_poly(lambda_1, x)])
    H_lambda = M.determinant()/(lambda_1+1-lambda_2)
    H_lambda = H_lambda.expand()
    H_lambda = polynomial(H_lambda, base_ring=QQ)
    return H_lambda


to_check = list(range(1, 1001))
to_check.extend([2401, 4375])
for n in to_check:
  if n%10 == 0:
    print(n)
  p = Partition([n, n])
  R_lambda = compute_wronskian_hermite(p)
  s = core_size(p)
  R_lambda = R_lambda / x^s
  R_lambda = R(R_lambda)
  if not R_lambda.is_irreducible():
    print(p, s, R_lambda)
