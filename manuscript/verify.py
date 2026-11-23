import itertools
import json
from dataclasses import dataclass
from fractions import Fraction
from math import factorial, gcd, prod


@dataclass
class Group:
    name: str
    table: tuple

    def __post_init__(self):
        self.order = len(self.table)
        self.elements = frozenset(range(self.order))
        assert all(self.table[0][element] == element == self.table[element][0]
                   for element in self.elements)
        self.inverses = tuple(next(inverse for inverse in self.elements
                                   if self.table[element][inverse] == 0
                                   and self.table[inverse][element] == 0)
                              for element in self.elements)
        assert all(self.table[self.table[first][second]][third]
                   == self.table[first][self.table[second][third]]
                   for first in self.elements for second in self.elements
                   for third in self.elements)

    def closure(self, generators):
        members = {0, *generators}
        while True:
            expanded = members | {self.table[first][second]
                                  for first in members for second in members}
            if expanded == members:
                return frozenset(members)
            members = expanded

    def subgroups(self):
        known = {frozenset({0})}
        pending = [frozenset({0})]
        while pending:
            subgroup = pending.pop()
            for element in sorted(self.elements - subgroup):
                generated = self.closure(subgroup | {element})
                if generated not in known:
                    known.add(generated)
                    pending.append(generated)
        return sorted(known, key=lambda subgroup: (len(subgroup), tuple(sorted(subgroup))))

    def normal(self, subgroup):
        return all(self.table[self.table[element][member]][self.inverses[element]]
                   in subgroup for element in self.elements for member in subgroup)

    def cosets(self, subgroup, left=True):
        remaining = set(self.elements)
        result = []
        while remaining:
            representative = min(remaining)
            coset = frozenset(self.table[representative][member] if left
                              else self.table[member][representative]
                              for member in subgroup)
            assert coset <= remaining
            result.append(coset)
            remaining -= coset
        return result

    def components(self, subgroup):
        remaining = set(self.elements)
        result = []
        while remaining:
            representative = min(remaining)
            component = frozenset(self.table[self.table[first][representative]][second]
                                  for first in subgroup for second in subgroup)
            assert component <= remaining
            assert len(component) % len(subgroup) == 0
            result.append((component, len(component) // len(subgroup)))
            remaining -= component
        return result

    def probability(self, subgroup):
        return prod((weight(size) for _, size in self.components(subgroup)), start=Fraction(1))

    def count_transversals(self, subgroup):
        left_cosets = self.cosets(subgroup)
        right_cosets = self.cosets(subgroup, left=False)
        right_index = {element: index for index, coset in enumerate(right_cosets)
                       for element in coset}
        choices = len(subgroup) ** len(left_cosets)
        assert choices <= 100000
        accepted = sum(len({right_index[element] for element in selected}) == len(left_cosets)
                       for selected in itertools.product(*(sorted(coset) for coset in left_cosets)))
        return Fraction(accepted, choices), choices

    def quotient(self, normal):
        assert self.normal(normal)
        cosets = self.cosets(normal)
        projection = {element: index for index, coset in enumerate(cosets) for element in coset}
        representatives = [min(coset) for coset in cosets]
        table = tuple(tuple(projection[self.table[first][second]] for second in representatives)
                      for first in representatives)
        return Group(self.name + '/N', table), projection


def weight(size):
    assert size >= 1
    return Fraction(factorial(size), size ** size)


def indexed_group(name, labels, multiply):
    labels = tuple(labels)
    indices = {label: index for index, label in enumerate(labels)}
    return Group(name, tuple(tuple(indices[multiply(first, second)] for second in labels)
                             for first in labels))


def cyclic(order):
    return indexed_group(f'C{order}', range(order), lambda first, second: (first + second) % order)


def dihedral(rotation_order):
    return indexed_group(f'D{2 * rotation_order}',
                         itertools.product(range(rotation_order), range(2)),
                         lambda first, second: ((first[0] + (-1) ** first[1] * second[0])
                                                % rotation_order, (first[1] + second[1]) % 2))


def permutations(degree, alternating=False):
    labels = list(itertools.permutations(range(degree)))
    if alternating:
        labels = [permutation for permutation in labels
                  if sum(permutation[first] > permutation[second]
                         for first in range(degree) for second in range(first + 1, degree)) % 2 == 0]
    return indexed_group(('A' if alternating else 'S') + str(degree), labels,
                         lambda first, second: tuple(first[second[index]] for index in range(degree)))


def direct_product(first, second):
    return indexed_group(first.name + 'x' + second.name,
                         itertools.product(range(first.order), range(second.order)),
                         lambda left, right: (first.table[left[0]][right[0]],
                                               second.table[left[1]][right[1]]))


def quaternion():
    basis_products = (((1, 0), (1, 1), (1, 2), (1, 3)),
                      ((1, 1), (-1, 0), (1, 3), (-1, 2)),
                      ((1, 2), (-1, 3), (-1, 0), (1, 1)),
                      ((1, 3), (1, 2), (-1, 1), (-1, 0)))

    def multiply(first, second):
        sign, basis = basis_products[first[1]][second[1]]
        return first[0] * second[0] * sign, basis

    return indexed_group('Q8', itertools.product((1, -1), range(4)), multiply)


def check_group(group):
    subgroups = group.subgroups()
    normals = [subgroup for subgroup in subgroups if group.normal(subgroup)]
    probabilities = {subgroup: group.probability(subgroup) for subgroup in subgroups}
    components = {subgroup: group.components(subgroup) for subgroup in subgroups}
    enumerated_choices = 0
    for subgroup in subgroups:
        actual, choices = group.count_transversals(subgroup)
        assert actual == probabilities[subgroup], (group.name, subgroup, actual)
        assert (actual == 1) == group.normal(subgroup)
        enumerated_choices += choices
    pairs = strict = equalities = component_surjections = quotient_checks = 0
    for normal in normals:
        quotient, projection = group.quotient(normal)
        for subgroup in subgroups:
            extended = frozenset(group.table[first][second] for first in subgroup for second in normal)
            assert extended in probabilities
            assert probabilities[subgroup] <= probabilities[extended]
            fibers = [[] for _ in components[extended]]
            for component, size in components[subgroup]:
                matching = [index for index, (downstairs, _) in enumerate(components[extended])
                            if component <= downstairs]
                assert len(matching) == 1
                target_index = matching[0]
                downstairs, downstairs_size = components[extended][target_index]
                image_cosets = {coset for coset in group.cosets(extended) if coset & component}
                target_cosets = {coset for coset in group.cosets(extended) if coset <= downstairs}
                assert image_cosets == target_cosets
                assert size >= downstairs_size
                fibers[target_index].append(size)
                component_surjections += 1
            assert all(fibers)
            equality_criterion = all(all(size == 1 for size in fiber) if downstairs_size == 1
                                     else fiber == [downstairs_size]
                                     for fiber, (_, downstairs_size) in zip(fibers, components[extended]))
            assert equality_criterion == (probabilities[subgroup] == probabilities[extended])
            assert equality_criterion == (normal <= subgroup or group.normal(subgroup))
            image = frozenset(projection[element] for element in subgroup)
            assert probabilities[subgroup] <= quotient.probability(image)
            assert (probabilities[subgroup] == quotient.probability(image)) == equality_criterion
            if normal <= subgroup:
                image = frozenset(projection[element] for element in subgroup)
                assert quotient.probability(image) == probabilities[subgroup]
                quotient_checks += 1
            pairs += 1
            strict += probabilities[subgroup] < probabilities[extended]
            equalities += probabilities[subgroup] == probabilities[extended]
    return {'group': group.name, 'order': group.order, 'subgroups': len(subgroups),
            'normal_subgroups': len(normals), 'subgroup_normal_pairs': pairs,
            'strict_inequalities': strict, 'equalities': equalities,
            'component_surjections': component_surjections,
            'quotient_identities': quotient_checks, 'enumerated_transversals': enumerated_choices}


def check_product(first, second):
    product_group = direct_product(first, second)
    first_subgroups = first.subgroups()
    second_subgroups = second.subgroups()
    product_subgroups = product_group.subgroups()
    probabilities = {subgroup: product_group.probability(subgroup)
                     for subgroup in product_subgroups}
    rectangular = set()
    rectangular_checks = 0
    for first_subgroup in first_subgroups:
        for second_subgroup in second_subgroups:
            subgroup = frozenset(first_element * second.order + second_element
                                 for first_element in first_subgroup
                                 for second_element in second_subgroup)
            rectangular.add(subgroup)
            first_probability = first.probability(first_subgroup)
            second_probability = second.probability(second_subgroup)
            assert probabilities[subgroup] >= (first_probability ** (second.order // len(second_subgroup))
                                                * second_probability ** (first.order // len(first_subgroup)))
            rectangular_checks += 1
    coprime = gcd(first.order, second.order) == 1
    if coprime:
        assert rectangular == set(product_subgroups)
        first_minimum = min(first.probability(subgroup) for subgroup in first_subgroups)
        second_minimum = min(second.probability(subgroup) for subgroup in second_subgroups)
        assert min(probabilities.values()) == min(first_minimum ** second.order,
                                                  second_minimum ** first.order)
    return {'group': product_group.name, 'order': product_group.order,
            'subgroups': len(product_subgroups), 'rectangular_checks': rectangular_checks,
            'coprime_minimum_identity': coprime,
            'minimum_probability': str(min(probabilities.values()))}


def main():
    groups = [cyclic(order) for order in (1, 2, 3, 4, 5, 6, 8, 12, 24)]
    groups += [dihedral(order) for order in (3, 4, 5, 6, 8, 10, 12)]
    groups += [quaternion(), permutations(4, alternating=True), permutations(4)]
    groups += [direct_product(cyclic(2), cyclic(2)),
               direct_product(direct_product(cyclic(2), cyclic(2)), cyclic(2)),
               direct_product(cyclic(3), cyclic(3)), direct_product(cyclic(2), cyclic(4)),
               direct_product(cyclic(2), cyclic(6))]
    rows = [check_group(group) for group in groups]
    symmetric = permutations(3)
    involution = next(element for element in range(1, symmetric.order)
                      if symmetric.table[element][element] == 0)
    middle = symmetric.closure({involution})
    assert [symmetric.probability(subgroup) for subgroup in
            [frozenset({0}), middle, symmetric.elements]] == [Fraction(1), Fraction(1, 2), Fraction(1)]
    for size in range(1, 100):
        assert weight(size + 1) / weight(size) == Fraction(size, size + 1) ** size
        assert 0 < weight(size + 1) < weight(size) <= 1
    for first_size in range(1, 25):
        for second_size in range(1, 25):
            assert weight(first_size * second_size) >= (weight(first_size) ** second_size
                                                        * weight(second_size) ** first_size)
    product_factors = [cyclic(1), cyclic(2), cyclic(3), cyclic(5), permutations(3), dihedral(4)]
    products = [check_product(first, second)
                for index, first in enumerate(product_factors)
                for second in product_factors[index:]
                if gcd(first.order, second.order) == 1 or first.order == second.order == 6]
    totals = {key: sum(row[key] for row in rows) for key in rows[0] if key not in ('group', 'order')}
    print(json.dumps({'status': 'PASS', 'groups': rows, 'totals': totals,
                      'weight_recurrence_range': [1, 99],
                      'grid_weight_range': {'first': [1, 24], 'second': [1, 24], 'pairs': 576},
                      'direct_products': products,
                      'structural_equality_checks': totals['subgroup_normal_pairs'],
                      'homomorphism_inequality_and_equality_checks': totals['subgroup_normal_pairs'],
                      'S3_subgroup_chain_probabilities': ['1', '1/2', '1'],
                      'scope': 'Exact finite checks of transversal counts, normal-extension inequalities, '
                               'double-coset and structural equality criteria, homomorphism consequences, '
                               'quotient identities, grid weights and direct-product bounds for the listed groups.'},
                     indent=2))


if __name__ == '__main__':
    main()
