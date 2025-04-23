# Room Placement Solver Using Graph Theory and MILP

This project solves the room placement problem in a square grid with three corners removed, maximizing the number of adjacent rooms that share walls. The problem is formulated using graph theory and solved using Mixed Integer Linear Programming (MILP) and constraint programming (CP-SAT) from Google's OR-Tools. The solution also visualizes the room layout using Matplotlib.

## Problem Overview
- The problem involves placing rooms (represented as nodes in a graph) inside a grid boundary while maximizing adjacencies between rooms. The layout must adhere to constraints such as room dimensions, grid boundary, and non-overlap between rooms.
- A mathematical model is formulated to represent room positions, dimensions, and adjacencies using variables for room placement and rotations.

## Features
- **Graph-Theoretic Interpretation**: The rooms are represented as physical rectangles, and edges (adjacencies) become shared walls.
- **MILP Model**: The solver uses OR-Tools' CP-SAT solver to maximize shared walls between adjacent rooms.
- **Visualization**: Matplotlib is used to visualize the final room layout.

## Installation

1. **Install Python**  
   Ensure that Python 3.7 or newer is installed. You can download Python from [here](https://www.python.org/downloads/).

2. **Install Required Packages**  
   After extracting the project files, open a terminal in the project directory and run:
   ```
   pip install -r requirements.txt
   ```
   This will install all necessary dependencies.

## Usage

1. To run the program, navigate to the project folder and execute:
   ```
   python gan_v2.py
   ```

2. Follow the on-screen menu to either:
   - Define your custom problem by specifying room sizes, adjacencies, and boundary.
   - Use built-in example scenarios.

The program will compute the room layout and display a visual representation of the solution.

## Key Components
- **Room Placement Algorithm**: Solves the optimization problem and visualizes the solution.
- **Graph-Theoretic Analysis**: The solution is analyzed for planarity and connectivity properties.
- **Constraints**: Various constraints such as room dimensions, non-overlap, and adjacency are modeled and enforced.

## Tips
- If you encounter issues with missing packages, make sure the requirements are installed as described above.
- Use `python3` instead of `python` if your system requires it.

## Conclusion
The algorithm effectively combines graph theory with constraint programming to solve a real-world spatial arrangement problem. By using Google's OR-Tools for optimization and Matplotlib for visualization, it provides an efficient solution for room placement and adjacency maximization in architectural and space allocation tasks.
