#! /usr/bin/env python
import getopt
import logging
import os
import pickle
import sys

from classifiers.pdfrate_online import get_fitness_func
from lib.common import deepcopy
from lib.common import setup_logging
from lib.common import visited_flag


# A good case:
# Classifier: pdfrate_bm, threshold=50
# Sandbox: gmu_bm, threshold=75

class GPPdf:
    def __init__(self,
                 job_dir,
                 seed_file_path,
                 logger,
                 ext_genome,
                 gp_params,
                 fitness_func,
                 sandbox_func,
                 hc_step
                 ):
        self.logger = logger
        self.job_dir = job_dir

        self.fitness_func = fitness_func
        self.sandbox_func = sandbox_func
        self.hc_step = hc_step

        # Load the seed.
        self.seed_file_path = seed_file_path
        self.seed_fitness = self.fitness([self.seed_file_path])[0]  # , self.seed_sha1
        self.seed_root = PdfGenome.load_genome(seed_file_path)
        self.logger.info("Loaded %s as PDF seed, fitness %.2f." % (seed_file_path, self.seed_fitness))

        # Load the external genome.
        self.ext_genome = ext_genome

        # Initiate some parameters.
        self.gp_params = gp_params
        self.pop_size = gp_params['pop_size']
        self.max_gen = gp_params['max_gen']
        self.mut_rate = gp_params['mut_rate']

    def save_trace_generation(self, variant, path):
        folder = "./variants/hc_%d/traces_generation" % self.hc_count
        folder = os.path.join(self.job_dir, folder)
        if not os.path.isdir(folder):
            os.makedirs(folder)
        path = os.path.join(folder, path)
        PdfGenome.save_to_file(variant, path)
        return path

    def save_variants_to_files(self):
        folder = "./variants/hc_%d/generation_%d" % (self.hc_count, self.generation)
        folder = os.path.join(self.job_dir, folder)
        if not os.path.isdir(folder):
            os.makedirs(folder)
        file_paths = []
        for j in range(len(self.popul)):
            path = "./variants/hc_%d/generation_%d/%d.pdf" % (self.hc_count, self.generation, j)
            path = os.path.join(self.job_dir, path)
            file_paths.append(path)
            PdfGenome.save_to_file(self.popul[j], path)
        return file_paths

    def fitness(self, *args):
        return self.fitness_func(*args)

    def run(self):
        self.logger.info("Start a gp task with %s" % (self.gp_params))

        score_file_name = os.path.join(self.job_dir, "fitness_scores.pickle")
        self.fitness_scores = {}
        self.hc_count = 0

        max_hc_count = 5

        while True:
            if self.hc_count > max_hc_count:
                self.logger.info("Giving up!")
                break

            self.generation = 0
            self.popul = []
            # Generate traces

            self.traces = []
            self.bin_search_range = []

            self.logger.info("HC step %d" % self.hc_count)

            generation_count = 0
            while len(self.traces) < self.pop_size:
                possible_new_mutation = PdfGenome.mutation(deepcopy(self.seed_root), self.mut_rate, self.ext_genome,
                                                           max_mut=2 ** self.max_gen)
                path = self.save_trace_generation(possible_new_mutation, "%d.pdf" % generation_count)
                score = self.fitness([path])[0]
                if not score:
                    self.popul.append(possible_new_mutation)
                    self.traces.append(possible_new_mutation.active_trace)
                    self.bin_search_range.append([0, len(possible_new_mutation.active_trace)])
                    test_pdf = Trace.generate_variant_from_trace(deepcopy(self.seed_root),
                                                                 possible_new_mutation.active_trace,
                                                                 self.ext_genome)
                    self.save_trace_generation(test_pdf, "%d_test.pdf" % generation_count)
                generation_count += 1
            self.save_variants_to_files()

            while self.generation < self.max_gen:
                self.generation += 1
                # Generate samples
                ends = []
                for i in range(self.pop_size):
                    end = (self.bin_search_range[i][0] + self.bin_search_range[i][1]) / 2
                    self.popul[i] = Trace.generate_variant_from_trace(deepcopy(self.seed_root), self.traces[i][:end],
                                                                      self.ext_genome)
                    ends.append(end)

                file_paths = self.save_variants_to_files()
                scores = self.fitness(file_paths)
                self.fitness_scores[self.generation] = scores
                pickle.dump(self.fitness_scores, open(score_file_name, 'wb'))
                self.logger.info("Fitness scores at generation %d: %s" % (self.generation, scores))

                for i in range(self.pop_size):
                    if scores[i]:
                        self.bin_search_range[i][0] = ends[i]
                    else:
                        self.bin_search_range[i][1] = ends[i]

            self.generation += 1
            for i in range(self.pop_size):
                end = self.bin_search_range[i][1]
                self.popul[i] = Trace.generate_variant_from_trace(deepcopy(self.seed_root), self.traces[i][:end],
                                                                  self.ext_genome)
            file_paths = self.save_variants_to_files()
            scores = self.fitness(file_paths)
            sandbox_result = self.sandbox_func(file_paths)
            flag = False
            for i in range(self.pop_size):
                if sandbox_result[i] and not scores[i]:
                    self.logger.info("Variant %d at gen %d passed!" % (i, self.generation))
                    flag = True
                    break
            if flag:
                break

            # nothing passed, bin search on sandbox
            classifier_ends = [self.bin_search_range[i][1] for i in range(self.pop_size)]

            self.bin_search_range = [[0, len(self.traces[i])] for i in range(self.pop_size)]
            while self.generation < 2 * self.max_gen:
                self.generation += 1
                # Generate samples
                ends = []
                for i in range(self.pop_size):
                    end = (self.bin_search_range[i][0] + self.bin_search_range[i][1]) / 2
                    self.popul[i] = Trace.generate_variant_from_trace(deepcopy(self.seed_root), self.traces[i][:end],
                                                                      self.ext_genome)
                    ends.append(end)

                file_paths = self.save_variants_to_files()
                scores = self.sandbox_func(file_paths)
                self.fitness_scores[self.generation] = scores
                pickle.dump(self.fitness_scores, open(score_file_name, 'wb'))
                self.logger.info("Sandbox at generation %d: %s" % (self.generation, scores))

                for i in range(self.pop_size):
                    if scores[i]:
                        self.bin_search_range[i][0] = ends[i]
                    else:
                        self.bin_search_range[i][1] = ends[i]

            sandbox_ends = [self.bin_search_range[i][0] for i in range(self.pop_size)]

            min_i = -1
            min_distance = 99999999
            for i in range(self.pop_size):
                this_distance = classifier_ends[i] - sandbox_ends[i]
                self.logger.info("Difference of %d: %d" % (i, this_distance))
                if this_distance < min_distance:
                    min_distance = this_distance
                    min_i = i
                if this_distance < 0:
                    self.logger.critical("========== Take a look at %d!" % i)
            self.logger.info("Picked variant %d at generation %d" % (min_i, self.generation))

            self.hc_count += 1
            new_root = Trace.generate_variant_from_trace(deepcopy(self.seed_root),
                                                         self.traces[min_i][:int(sandbox_ends[min_i] * hc_step)],
                                                         self.ext_genome)
            new_root_path = self.save_trace_generation(new_root, 'root.pdf')
            self.seed_root = PdfGenome.load_genome(new_root_path)
        return True


def get_opt(argv):
    start_file = None
    ext_genome_folder = None
    pop_size = None
    max_gen = None
    mut_rate = None
    round_id = '1'
    field = 'contagio_bm'
    threshold = 50
    hc_step = 0.75

    # try:
    opts, args = getopt.getopt(argv[1:], "s:e:p:g:m:r:f:t:h:", [  # "classifier=",
        "sfile=",
        "extgenome=",
        "popu=",
        "gen=",
        "mut=",
        "rid=",
        "field=",
        "threshold=",
        "hc_step="
    ])

    for opt, arg in opts:
        # elif opt in ("-c", "--classifier"):
        #     classifier_name = arg
        if opt in ("-s", "--sfile"):
            start_file = arg
        elif opt in ("-e", "--extgenome"):
            ext_genome_folder = arg
        elif opt in ("-p", "--popu"):
            pop_size = int(arg)
        elif opt in ("-g", "--gen"):
            max_gen = int(arg)
        elif opt in ("-m", "--mut"):
            mut_rate = float(arg)
        elif opt in ("-r", "--rid"):
            round_id = arg
        elif opt in ("-f", "--field"):
            field = arg
        elif opt in ("-t", "--threshold"):
            threshold = float(arg)
        elif opt in ("-h", "--hc_step"):
            hc_step = float(arg)

    return start_file, ext_genome_folder, pop_size, max_gen, mut_rate, round_id, field, threshold, hc_step


if __name__ == "__main__":
    start_file_path, ext_genome_folder, pop_size, max_gen, mut_rate, round_id, field, threshold, hc_step = get_opt(
        sys.argv)

    job_dir = "./results/log_r%s" % round_id
    if not os.path.isdir(job_dir):
        os.makedirs(job_dir)

    log_file_path = os.path.join(job_dir, visited_flag)
    setup_logging(log_file_path)
    logger = logging.getLogger('gp.core')
    logger.info("Starting logging for a GP process...")

    # Due to logging is called in pdfrw, they have to be imported after basicConfig of logging.
    # Otherwise, the above basicConfig would be overridden.
    from lib.pdf_genome import PdfGenome
    from lib.trace import Trace

    # if classifier_name == 'pdfrate':
    #     from lib.fitness import fitness_pdfrate as fitness_func
    # elif classifier_name == 'hidost':
    #     from lib.fitness import fitness_hidost as fitness_func
    # elif classifier_name == "hidost_pdfrate":
    #     from lib.fitness import fitness_hidost_pdfrate as fitness_func
    # elif classifier_name == "hidost_pdfrate_mean":
    #     from lib.fitness import fitness_hidost_pdfrate_mean as fitness_func
    # elif classifier_name == "hidost_pdfrate_sigmoid":
    #     from lib.fitness import fitness_hidost_pdfrate_sigmoid as fitness_func

    gp_params = {'pop_size': pop_size, 'max_gen': max_gen, 'mut_rate': mut_rate}
    ext_genome = PdfGenome.load_external_genome(ext_genome_folder)

    fitness_func = get_fitness_func(field, threshold)

    gp = GPPdf(job_dir=job_dir,
               seed_file_path=start_file_path,
               logger=logger,
               ext_genome=ext_genome,
               gp_params=gp_params,
               fitness_func=fitness_func,
               sandbox_func=get_fitness_func('gmu_bm', threshold=65),
               hc_step=hc_step
               )
    gp.run()
