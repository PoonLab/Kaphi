# This file is part of Kaphi.

# Kaphi is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Kaphi is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Kaphi.  If not, see <http://www.gnu.org/licenses/>.

require(Kaphi, quietly=TRUE)
require(RUnit, quietly=TRUE)
source('tests/fixtures/simple-trees.R')   # not sure if this is needed

test.compartmental.model <- function() {
  checkException(result <- )
  
  result <-
  expected <-
  checkEquals(expected, result)
}