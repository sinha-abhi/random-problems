/*
 * Fully sorting a k-storted array.
 *
 * Abhi Sinha, sinha45@purdue.edu
 * January 28, 2020
 */

#include <iostream>
#include <queue>
#include <vector>

std::vector<int> ksort(std::vector<int> &array, int k) {
    std::vector<int> sorted;
    std::priority_queue<int, std::vector<int>, std::greater<int>> q;

    for (int i = 0; i < k + 1; i++)
        q.push(array[i]);

    for (int i = k + 1; i < array.size(); i++) {
        sorted.push_back(q.top());
        q.pop();
        q.push(array[i]);
    }

    while (!q.empty()) {
        sorted.push_back(q.top());
        q.pop();
    }

    return sorted;
}

void ksort_modify(std::vector<int> &array, int k) {
    std::priority_queue<int, std::vector<int>, std::greater<int>> q;

    for (int i = 0; i < k + 1; i++)
        q.push(array[i]);

    int pos = 0;
    for (int i = k + 1; i < array.size(); i++) {
        array[pos++] = q.top();
        q.pop();
        q.push(array[i]);
    }

    while (!q.empty()) {
        array[pos++] = q.top();
        q.pop();
    }
}
