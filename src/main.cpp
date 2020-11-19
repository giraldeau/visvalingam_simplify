//
//
// 2013 (c) Mathieu Courtemanche
//
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include "visvalingam_algorithm.h"
#include "ogr.h"
#include "csv.h"
#include "heap.hpp"
#include "decimate.h"

int export_filtered_csv(const std::string &filename, const std::string &suffix, std::vector<int> idx)
{
    std::ofstream os_filtered(std::string(filename) + suffix);
    std::ifstream is(filename);
    if (!is.is_open()) {
        return 1;
    }

    int lineno = 0;
    std::string line;
    size_t curr = 0;
    while(std::getline(is, line)) {
        if (idx[curr] == lineno) {
            os_filtered << line << "\n";
            curr++;
        }
        lineno++;
    }
    return 0;
}

int process_csv(const char *filename, bool print_source, size_t ratio, const std::vector<int> &cols, int &keep_col)
{
    Linestring shape;
    Linestring shape_simplified;
    Linestring shape_decimated;
    std::vector<bool> keep_nodes;
    std::ifstream is(filename);

    if (!is.is_open()) {
        return 1;
    }

    std::cout << "using columns: ";
    std::copy(cols.begin(), cols.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;

    CSVIterator iter(is);
    int lineno = 1;
    int keepcount = 0;

    while(iter != CSVIterator()) {
       std::vector<double> coords(3, 0.0);
       bool keep = false;

       const CSVRow &row = *iter;
        for (size_t col = 0; col < cols.size(); ++col) {
            size_t col_id = cols[col];
            if (col_id < row.size()) {
                coords[col] = std::stod(row[col_id]);
            } else {
                std::cout << "error: failed to parse column " << col_id << " at line " << lineno << std::endl;
                return -1;
            }
        }
        if (keep_col >= 0 && keep_col < row.size()) {
            keep = std::stoi(row[keep_col]) == 0 ? false : true;
            if (keep) {
                keepcount++;
            }
        }

        shape.push_back(Point(coords[0], coords[1], coords[2]));
        keep_nodes.push_back(keep);

        iter++;
        lineno++;
    }
    std::cout << "number of points: " << shape.size() << std::endl;
    std::cout << "number of lines : " << lineno - 1 << std::endl;
    std::cout << "number of keep  : " << keepcount << std::endl;

    VertexFilter keepFilter = [&keep_nodes](const VertexIndex &i) -> bool { return keep_nodes[i]; };
    Visvalingam_Algorithm vis_algo(shape);
    double threshold = vis_algo.area_threshold_for_ratio(ratio);
    std::vector<int> idx;
    std::vector<int> idx2;
    vis_algo.simplify(threshold, &shape_simplified, &idx, keepFilter);

    decimate(shape, ratio, shape_decimated, idx2);

    std::cout << "area threshold:   " << threshold << std::endl;
    std::cout << "original shape:   " << shape.size() << " points" << std::endl;
    std::cout << "simplified shape: " << shape_simplified.size() << " points" << std::endl;
    std::cout << "decimateed shape: " << shape_decimated.size() << " points" << std::endl;


    std::ofstream os_pts(std::string(filename) + ".out");
    for (const Point &p : shape_simplified) {
        os_pts << p.X << "," << p.Y << "," << p.Z << "\n";
    }

    std::ofstream os_areas(std::string(filename) + ".areas");
    for (int i = 1; i < shape.size() - 1; i++) {
        os_areas << effective_area(shape[i], shape[i-1], shape[i+1]) << "\n";
    }

    export_filtered_csv(filename, ".filtered", idx);
    export_filtered_csv(filename, ".decimated", idx);

    return 0;
}

enum InputFormat {
    FORMAT_OGR,
    FORMAT_CSV,
};

int main(int argc, char **argv)
{
    bool print_source = false;
    const char* filename = NULL;
    size_t ratio = 50;
    int filter_col = -1;
    std::vector<int> coord_cols{0, 1, 2};
    InputFormat file_format = FORMAT_OGR;
    int res = 0;
    for (int i=1; i < argc; ++i)
    {
        if (strcmp(argv[i], "--file") == 0 && (i+1) < argc)
        {
            ++i;
            filename = argv[i];
        }
        else if (strcmp(argv[i], "--ratio") == 0 && (i+1) < argc)
        {
            ++i;
            ratio = static_cast<size_t>(std::atoi(argv[i]));
            std::cout << "using ratio: " << ratio << std::endl;
        }
        else if (strcmp(argv[i], "--filter-col") == 0 && (i+1) < argc)
        {
            ++i;
            filter_col = static_cast<size_t>(std::atoi(argv[i]));
        }
        else if (strcmp(argv[i], "--cols") == 0 && (i+1) < argc)
        {
            ++i;
            coord_cols.clear();
            std::string cols_arg(argv[i]);
            std::string item;
            std::stringstream ss(cols_arg);
            while (std::getline(ss, item, ',')) {
                coord_cols.push_back(std::stoi(item));
            }
        }
        else if (strcmp(argv[i], "--dump-source") == 0)
        {
            print_source = true;
        }
        else if (strcmp(argv[i], "--format=ogr") == 0)
        {
            file_format = FORMAT_OGR;
        }
        else if (strcmp(argv[i], "--format=csv") == 0)
        {
            file_format = FORMAT_CSV;
        }
        else
        {
            std::cout << "unkown argument: " << argv[i] << std::endl;
            return 1;
        }
    }

    if (!filename)
    {
        std::cout << "error: input file required" << std::endl;
        return 1;
    }

    if (ratio <= 0 || ratio >= 100)
    {
        std::cout << "error: bad ratio, must be between 1 and 99, got: " << ratio << std::endl;
        return 1;
    }

    switch(file_format)
    {
    case FORMAT_OGR:
        std::cout << "processing_ogr" << std::endl;
        res = process_ogr(filename, print_source);
        break;
    case FORMAT_CSV:
        std::cout << "processing_csv" << std::endl;
        res = process_csv(filename, print_source, ratio, coord_cols, filter_col);
        break;
    }

    if (res != 0)
    {
        std::cout << "error: input file loading failed" << std::endl;
        return res;
    }

    return res;
}
