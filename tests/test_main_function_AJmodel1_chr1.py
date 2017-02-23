import json
import unittest

from main_function_AJmodel1_chr1 import main

class MyMain(unittest.TestCase):
    def test_basic(self):
        arguments  = ['',1,'ill_650_test.bed',1000,1278,'prior',1]
        results = main(arguments)
        print results

        def dump_results(results_to_dump):
            """Dump the results to a JSON file for comparison"""
            import json
            import time
            with open('tests/test_data/my_main_test_basic_result_{}.json'.format(time.time()), 'w') as fp:
                json.dump(results_to_dump, fp, indent=2, ensure_ascii=True)

        #dump_results(results)

        with open('tests/test_data/my_main_test_basic_result_1487808470.62.json', 'r') as fp:
            expected_results = json.load(fp, encoding='ascii')

        self.maxDiff = None
        self.assertListEqual(results, expected_results)


if __name__ == '__main__':
    unittest.main()
