/**
 * @file
 * @brief Implementation of the [Smith-Waterman
 * algorithm](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
 * for local sequence alignment
 *
 * @details
 * The Smith-Waterman algorithm is a dynamic programming algorithm for
 * determining similar regions between two sequences (nucleotide or protein
 * sequences). It performs local sequence alignment, which identifies the best
 * matching subsequence rather than aligning entire sequences.
 *
 * The algorithm works by:
 * 1. Creating a scoring matrix where each cell represents the maximum
 *    alignment score ending at that position
 * 2. Using match, mismatch, and gap penalties to calculate scores
 * 3. Allowing scores to reset to 0 (ensuring local rather than global
 *    alignment)
 * 4. Tracing back from the highest scoring position to reconstruct the
 *    alignment
 *
 * ### Algorithm
 * Given two sequences S1 and S2 of lengths m and n:
 * - Initialize a (m+1) × (n+1) matrix H with H[i][0] = H[0][j] = 0
 * - For each cell H[i][j]:
 *   - H[i][j] = max(0, H[i-1][j-1] + match/mismatch_score,
 *                   H[i-1][j] + gap_penalty, H[i][j-1] + gap_penalty)
 * - Find the maximum value in H matrix
 * - Traceback from maximum value to cell with value 0
 *
 * Time Complexity: O(m × n) where m and n are the lengths of the sequences
 * Space Complexity: O(m × n) for the scoring matrix
 *
 * @author [Ali Alimohammadi](https://github.com/AliAlimohammadi)
 */

#include <algorithm>   /// for std::max
#include <cassert>     /// for assert
#include <iostream>    /// for IO operations
#include <string>      /// for std::string
#include <vector>      /// for std::vector

/**
 * @namespace dynamic_programming
 * @brief Dynamic Programming algorithms
 */
namespace dynamic_programming {
/**
 * @namespace smith_waterman
 * @brief Functions for the Smith-Waterman algorithm
 */
namespace smith_waterman {

/**
 * @brief Calculate the score for a character pair
 *
 * @param a First character
 * @param b Second character
 * @param match_score Score for matching characters (typically positive)
 * @param mismatch_score Score for mismatching characters (typically negative)
 * @param gap_score Penalty for gaps (typically negative)
 * @return int The calculated score
 */
int score_function(char a, char b, int match_score, int mismatch_score,
                   int gap_score) {
    if (a == '-' || b == '-') {
        return gap_score;
    } else if (a == b) {
        return match_score;
    } else {
        return mismatch_score;
    }
}

/**
 * @brief Perform Smith-Waterman local sequence alignment
 *
 * @param query First sequence
 * @param subject Second sequence
 * @param match_score Score for matching characters (default: 1)
 * @param mismatch_score Score for mismatching characters (default: -1)
 * @param gap_score Penalty for gaps (default: -2)
 * @return std::vector<std::vector<int>> The scoring matrix
 */
std::vector<std::vector<int>> smith_waterman(const std::string& query,
                                               const std::string& subject,
                                               int match_score = 1,
                                               int mismatch_score = -1,
                                               int gap_score = -2) {
    int m = query.length();
    int n = subject.length();

    // Initialize scoring matrix with zeros
    std::vector<std::vector<int>> score(m + 1, std::vector<int>(n + 1, 0));

    // Fill the scoring matrix using dynamic programming
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            // Calculate score for match/mismatch
            int match_or_mismatch = score[i - 1][j - 1] +
                                     score_function(query[i - 1], subject[j - 1],
                                                    match_score, mismatch_score,
                                                    gap_score);

            // Calculate score for deletion (gap in subject)
            int delete_score = score[i - 1][j] + gap_score;

            // Calculate score for insertion (gap in query)
            int insert_score = score[i][j - 1] + gap_score;

            // Take maximum of all options, but never go below 0 (local
            // alignment)
            score[i][j] =
                std::max({0, match_or_mismatch, delete_score, insert_score});
        }
    }

    return score;
}

/**
 * @brief Perform traceback to reconstruct the optimal alignment
 *
 * @param score The score matrix from smith_waterman function
 * @param query Original query sequence
 * @param subject Original subject sequence
 * @param match_score Score for matching characters (default: 1)
 * @param mismatch_score Score for mismatching characters (default: -1)
 * @param gap_score Penalty for gaps (default: -2)
 * @return std::pair<std::string, std::string> The aligned sequences
 */
std::pair<std::string, std::string> traceback(
    const std::vector<std::vector<int>>& score, const std::string& query,
    const std::string& subject, int match_score = 1, int mismatch_score = -1,
    int gap_score = -2) {
    // Find the cell with maximum score
    int max_value = 0;
    int i_max = 0;
    int j_max = 0;

    for (size_t i = 0; i < score.size(); ++i) {
        for (size_t j = 0; j < score[i].size(); ++j) {
            if (score[i][j] > max_value) {
                max_value = score[i][j];
                i_max = i;
                j_max = j;
            }
        }
    }

    // If no significant alignment found, return empty strings
    if (i_max == 0 || j_max == 0) {
        return {"", ""};
    }

    // Traceback from the maximum scoring cell
    std::string align1;
    std::string align2;
    int i = i_max;
    int j = j_max;

    // Continue tracing back until we hit a cell with score 0
    while (i > 0 && j > 0 && score[i][j] > 0) {
        int current_score = score[i][j];

        // Check if we came from diagonal (match/mismatch)
        if (current_score ==
            score[i - 1][j - 1] +
                score_function(query[i - 1], subject[j - 1], match_score,
                               mismatch_score, gap_score)) {
            align1 = query[i - 1] + align1;
            align2 = subject[j - 1] + align2;
            --i;
            --j;
        }
        // Check if we came from above (deletion/gap in subject)
        else if (current_score == score[i - 1][j] + gap_score) {
            align1 = query[i - 1] + align1;
            align2 = '-' + align2;
            --i;
        }
        // Otherwise we came from left (insertion/gap in query)
        else {
            align1 = '-' + align1;
            align2 = subject[j - 1] + align2;
            --j;
        }
    }

    return {align1, align2};
}

}  // namespace smith_waterman
}  // namespace dynamic_programming

/**
 * @brief Self-test implementations
 * @returns void
 */
static void test() {
    // Test 1: Simple exact match
    auto score1 = dynamic_programming::smith_waterman::smith_waterman("AGT", "AGT");
    assert(score1[3][3] == 3);  // Perfect match should score 3
    std::cout << "Test 1 passed: Simple exact match\n";

    // Test 2: Partial match
    auto score2 = dynamic_programming::smith_waterman::smith_waterman("ACAC", "CA");
    assert(score2[3][2] == 2);  // Best local alignment scores 2
    std::cout << "Test 2 passed: Partial match\n";

    // Test 3: Traceback test
    auto result3 = dynamic_programming::smith_waterman::traceback(score2, "ACAC", "CA");
    assert(result3.first == "CA");
    assert(result3.second == "CA");
    std::cout << "Test 3 passed: Traceback\n";

    // Test 4: No match
    auto score4 = dynamic_programming::smith_waterman::smith_waterman("AAA", "TTT");
    int max_score = 0;
    for (const auto& row : score4) {
        max_score = std::max(max_score, *std::max_element(row.begin(), row.end()));
    }
    assert(max_score == 0);  // No matches should score 0
    std::cout << "Test 4 passed: No match\n";

    // Test 5: Empty strings
    auto score5 = dynamic_programming::smith_waterman::smith_waterman("", "AGT");
    assert(score5.size() == 1);
    assert(score5[0].size() == 4);
    std::cout << "Test 5 passed: Empty query\n";

    // Test 6: Longer sequences with gaps
    auto score6 = dynamic_programming::smith_waterman::smith_waterman("AGCT", "AGT");
    auto result6 = dynamic_programming::smith_waterman::traceback(score6, "AGCT", "AGT");
    // The alignment should handle the gap appropriately
    assert(!result6.first.empty());
    assert(!result6.second.empty());
    std::cout << "Test 6 passed: Longer sequences with gaps\n";

    // Test 7: Custom scoring
    auto score7 = dynamic_programming::smith_waterman::smith_waterman(
        "ACGT", "ACGT", 2, -1, -1);  // Higher match score
    assert(score7[4][4] == 8);  // 4 matches × 2 = 8
    std::cout << "Test 7 passed: Custom scoring\n";

    // Test 8: Case sensitivity (algorithm is case-sensitive)
    auto score8 = dynamic_programming::smith_waterman::smith_waterman("ACT", "act");
    int max_score8 = 0;
    for (const auto& row : score8) {
        max_score8 = std::max(max_score8, *std::max_element(row.begin(), row.end()));
    }
    assert(max_score8 == 0);  // Different cases, no match
    std::cout << "Test 8 passed: Case sensitivity\n";

    std::cout << "All tests passed successfully!\n";
}

/**
 * @brief Main function - only enabled for standalone testing
 * @returns 0 on exit
 */
int main() {
    test();  // Run self-test implementations

    // Example usage
    std::cout << "\n--- Example Usage ---\n";
    std::string query = "ACACACTA";
    std::string subject = "AGCACACA";

    std::cout << "Query:   " << query << "\n";
    std::cout << "Subject: " << subject << "\n\n";

    // Perform alignment
    auto score_matrix =
        dynamic_programming::smith_waterman::smith_waterman(query, subject);
    auto alignment =
        dynamic_programming::smith_waterman::traceback(score_matrix, query, subject);

    std::cout << "Aligned sequences:\n";
    std::cout << "Query:   " << alignment.first << "\n";
    std::cout << "Subject: " << alignment.second << "\n";

    // Find the maximum score
    int max_score = 0;
    for (const auto& row : score_matrix) {
        max_score = std::max(max_score, *std::max_element(row.begin(), row.end()));
    }
    std::cout << "Alignment score: " << max_score << "\n";

    return 0;
}
