#include "task_manager.hxx"

#include <cassert>

namespace CarpetX {

task_manager::task_manager() {}

task_manager::~task_manager() { assert(tasks.empty()); }

void task_manager::submit_serially(task_t task) {
  tasks.push_back(std::move(task));
}

void task_manager::submit(task_t task) {
#pragma omp critical(CarpetX_task_manager_submit)
  tasks.push_back(std::move(task));
}

void task_manager::run_tasks_serially() {
  for (const auto &task : tasks)
    task();
  tasks.clear();
}

void task_manager::run_tasks() {
  const std::size_t ntasks = tasks.size();
  if (ntasks == 0)
    return;
#pragma omp parallel for schedule(dynamic)
  for (std::size_t n = 0; n < ntasks; ++n)
    tasks[n]();
  tasks.clear();
}

} // namespace CarpetX
