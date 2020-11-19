#include <csv.h>
#include <polyline.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCylinderSource.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkNamedColors.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include <array>
#include <sstream>
#include <vector>

struct args {
  std::string filename;
  size_t ratio = 50;
  int filter_col = -1;
  std::vector<int> coord_cols{0, 1, 2};
};

int parse_args(int argc, char *argv[], struct args &args) {
  int res = 0;
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--file") == 0 && (i + 1) < argc) {
      ++i;
      args.filename = argv[i];
    } else if (strcmp(argv[i], "--ratio") == 0 && (i + 1) < argc) {
      ++i;
      args.ratio = static_cast<size_t>(std::atoi(argv[i]));
    } else if (strcmp(argv[i], "--filter-col") == 0 && (i + 1) < argc) {
      ++i;
      args.filter_col = static_cast<size_t>(std::atoi(argv[i]));
    } else if (strcmp(argv[i], "--cols") == 0 && (i + 1) < argc) {
      ++i;
      args.coord_cols.clear();
      std::string cols_arg(argv[i]);
      std::string item;
      std::stringstream ss(cols_arg);
      while (std::getline(ss, item, ',')) {
        args.coord_cols.push_back(std::stoi(item));
      }
    } else {
      std::cout << "unkown argument: " << argv[i] << std::endl;
      res = -1;
    }
  }
  return res;
}

/*
 * FIXME: de-duplicate the code with simplify
 */
int parse_csv(struct args &args, Linestring &shape, std::vector<bool> &keep_nodes) {
  std::ifstream ifs(args.filename);
  std::cout << "opening file: " << args.filename << std::endl;
  std::cout << ifs.is_open() << std::endl;
  if (!ifs.is_open()) {
    std::cout << "failed to open file " << args.filename;
    return -1;
  }

  std::cout << "using columns: ";
  std::copy(args.coord_cols.begin(), args.coord_cols.end(),
            std::ostream_iterator<int>(std::cout, " "));
  std::cout << std::endl;

  CSVIterator iter(ifs);
  int lineno = 1;
  int keepcount = 0;

  while (iter != CSVIterator()) {
    std::vector<double> coords(3, 0.0);
    bool keep = false;

    const CSVRow &row = *iter;
    for (size_t col = 0; col < args.coord_cols.size(); ++col) {
      size_t col_id = args.coord_cols[col];
      if (col_id < row.size()) {
        coords[col] = std::stod(row[col_id]);
      } else {
        std::cout << "error: failed to parse column " << col_id << " at line " << lineno
                  << std::endl;
        return -1;
      }
    }
    if (args.filter_col >= 0 && args.filter_col < row.size()) {
      keep = std::stoi(row[args.filter_col]) == 0 ? false : true;
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

  return 0;
}

void linestring_to_polydata(Linestring &shape, vtkPolyData *linesPolyData) {
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

  points->SetNumberOfPoints(shape.size());
  lines->SetNumberOfCells(shape.size() - 1);

  for (size_t i = 0; i < shape.size(); i++) {
    const Point &p = shape[i];
    points->SetPoint(i, p.X, p.Y, p.Z);
  }

  for (size_t i = 0; i < shape.size() - 1; i++) {
    lines->InsertNextCell(2);
    lines->InsertCellPoint(i);
    lines->InsertCellPoint(i + 1);
  }

  linesPolyData->SetPoints(points);
  linesPolyData->SetLines(lines);
}

int main(int argc, char *argv[]) {
  struct args args;

  if (parse_args(argc, argv, args) < 0) {
    std::cout << "error parsing arguments" << std::endl;
    return EXIT_FAILURE;
  }

  Linestring shape;
  std::vector<bool> keep_nodes;

  if (parse_csv(args, shape, keep_nodes) < 0) {
    std::cout << "error parsing csv file" << std::endl;
    return EXIT_FAILURE;
  }

  vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();
  linestring_to_polydata(shape, linesPolyData);

  // Set the background color.
  vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();
  std::array<unsigned char, 4> bkg{{255, 255, 255, 255}};
  colors->SetColor("BkgColor", bkg.data());

  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(linesPolyData);

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(colors->GetColor4d("Black").GetData());

  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->AddActor(actor);
  renderer->SetBackground(colors->GetColor3d("BkgColor").GetData());
  renderer->ResetCamera();
  renderer->GetActiveCamera()->Zoom(1.5);

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->SetSize(1024, 768);
  renderWindow->AddRenderer(renderer);
  renderWindow->SetWindowName("Path Viewer");

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor
      = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style
      = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle(style);

  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
