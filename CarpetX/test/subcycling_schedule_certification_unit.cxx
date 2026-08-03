#define CARPETX_SUBCYCLING_SCHEDULE_CERTIFICATION_STANDALONE 1

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using CCTK_INT = int;
struct cGH {};
enum cLanguage { LangNone, LangC, LangFortran };
enum cFunctionType { FunctionNoArgs, FunctionOneArg, FunctionStandard };

struct RDWR_entry {
  int varindex;
  int timelevel;
  int where_wr;
  int where_rd;
  int where_inv;
};

struct cFunctionData {
  cLanguage language;
  int (*FortranCaller)(cGH *, void *);
  cFunctionType type;
  int n_SyncGroups;
  int *SyncGroups;
  int meta;
  int meta_early;
  int meta_late;
  int global;
  int global_early;
  int global_late;
  int level;
  int singlemap;
  int local;
  int loop_meta;
  int loop_global;
  int loop_level;
  int loop_singlemap;
  int loop_local;
  int tags;
  int n_TriggerGroups;
  int *TriggerGroups;
  int n_WritesClauses;
  const char **WritesClauses;
  int n_ReadsClauses;
  const char **ReadsClauses;
  int n_InvalidatesClauses;
  const char **InvalidatesClauses;
  int n_RDWR;
  RDWR_entry *RDWR;
  char *where;
  char *routine;
  char *thorn;
};

struct cScheduleInventoryScope {
  unsigned int abi_version;
  std::size_t struct_size;
  std::uint64_t item_ordinal;
  std::uint64_t parent_local_ordinal;
  const char *description;
  const char *implementation;
  const cFunctionData *function_data;
  int n_storage_groups;
  const int *storage_groups;
  const int *storage_timelevels;
  int n_communication_groups;
  const int *communication_groups;
  int n_ifs;
  const char *const *ifs;
  int n_whiles;
  const char *const *whiles;
};

struct cScheduleInventoryRecord {
  unsigned int abi_version;
  std::size_t struct_size;
  const char *root_group;
  int entry_exit_included;
  std::uint64_t traversal_ordinal;
  std::uint64_t item_ordinal;
  std::uint64_t parent_local_ordinal;
  const char *description;
  const char *implementation;
  const cFunctionData *function_data;
  int n_storage_groups;
  const int *storage_groups;
  const int *storage_timelevels;
  int n_communication_groups;
  const int *communication_groups;
  int n_ifs;
  const char *const *ifs;
  int n_whiles;
  const char *const *whiles;
  std::size_t n_enclosing_groups;
  const cScheduleInventoryScope *enclosing_groups;
};

using cScheduleInventoryCallback =
    int (*)(const cScheduleInventoryRecord *, void *, void *);

struct cGroup {
  int grouptype;
  int vartype;
  int disttype;
  int dim;
  int numvars;
  int numtimelevels;
  int vectorgroup;
  int vectorlength;
  int tagstable;
  int centeringtable;
};

constexpr unsigned int CCTK_SCHEDULE_INVENTORY_ABI_VERSION = 1U;
constexpr int CCTK_SCHEDULE_INVENTORY_GROUP_NOT_FOUND = -1;
constexpr int CCTK_VARIABLE_INT = 120;
constexpr int CCTK_VARIABLE_REAL = 130;
constexpr int CCTK_VARIABLE_STRING = 151;
constexpr int CCTK_VARIABLE_POINTER = 160;
constexpr int CCTK_VARIABLE_POINTER_TO_CONST = 161;
constexpr int CCTK_VARIABLE_FPOINTER = 162;
constexpr int CCTK_SCALAR = 401;
constexpr int CCTK_GF = 402;
constexpr int CCTK_ARRAY = 403;
constexpr int CCTK_VALID_GHOSTS = 1;
constexpr int CCTK_VALID_BOUNDARY = 2;
constexpr int CCTK_VALID_INTERIOR = 4;
constexpr int UTIL_TABLE_FLAGS_CASE_INSENSITIVE = 1;

extern "C" {
int CCTK_ScheduleInventoryDry(const char *, cScheduleInventoryCallback, void *);
int CCTK_CallFunction(void *, cFunctionData *, cGH *);
const char *CCTK_FullGroupName(int);
int CCTK_GroupIndex(const char *);
const char *CCTK_FullVarName(int);
int CCTK_VarIndex(const char *);
int CCTK_GroupIndexFromVarI(int);
int CCTK_GroupData(int, cGroup *);
int CCTK_VarTypeI(int);
int CCTK_DeclaredTimeLevelsVI(int);
int CCTK_VarTypeSize(int);
int Util_TableQueryFlags(int);
int Util_TableQueryNKeys(int);
int Util_TableQueryMaxKeyLength(int);
int Util_TableItCreate(int);
int Util_TableItDestroy(int);
int Util_TableItQueryKeyValueInfo(int, int, char *, CCTK_INT *, CCTK_INT *);
int Util_TableGetGenericArray(int, int, int, void *, const char *);
int Util_TableItAdvance(int);
int Util_TableItQueryIsNull(int);
int Util_TableClone(int);
int Util_TableDestroy(int);
}

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsubobject-linkage"
#endif
#include "../src/subcycling_schedule_certification.cxx"
#include "../src/subcycling_schedule_certification_native_gate.hxx"
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace fake {

struct TagEntry {
  std::string key;
  int type{};
  int count{};
  std::vector<std::uint8_t> bytes;
  std::optional<int> get_result;
};

struct Table {
  int flags{UTIL_TABLE_FLAGS_CASE_INSENSITIVE};
  std::vector<TagEntry> entries;
  bool clone{};
};

struct Iterator {
  int table{-1};
  std::size_t index{};
};

enum class TableFailure {
  none,
  iterator_create,
  iterator_query,
  iterator_advance,
  iterator_end_state,
  clone,
};

struct GroupInfo {
  std::string name;
  cGroup data{};
};

struct VarInfo {
  std::string name;
  int group{};
  int type{};
  int timelevels{};
};

static std::map<int, Table> tables;
static std::map<int, Iterator> iterators;
static int next_table_handle = 100;
static int next_iterator_handle = 1000;
static int live_clones = 0;
static int live_iterators = 0;
static int poison_invocations = 0;
static int &live_cloned_table_count = live_clones;
static int &poison_invocation_count = poison_invocations;
static bool throw_on_table_get = false;
static TableFailure table_failure = TableFailure::none;
static bool break_group_roundtrip = false;
static bool break_var_roundtrip = false;
static std::vector<std::string> inventory_calls;
static std::vector<GroupInfo> group_infos;
static std::vector<VarInfo> var_infos;

template <class T> std::vector<std::uint8_t> bytes_of(const T &value) {
  std::vector<std::uint8_t> bytes(sizeof(T));
  std::memcpy(bytes.data(), &value, sizeof(T));
  return bytes;
}

static int make_table(const bool reverse = false) {
  const int integer_value = 9;
  const double real_value = 1.5;
  Table table;
  table.entries = {{"Beta", CCTK_VARIABLE_INT, 1, bytes_of(integer_value), {}},
                   {"alpha", CCTK_VARIABLE_REAL, 1, bytes_of(real_value), {}}};
  if (reverse)
    std::reverse(table.entries.begin(), table.entries.end());
  const int handle = next_table_handle++;
  tables.emplace(handle, std::move(table));
  return handle;
}

struct FunctionFixture {
  cFunctionData data{};
  std::vector<int> sync_groups;
  std::vector<int> trigger_groups;
  std::vector<std::string> reads{"Test::u"};
  std::vector<std::string> writes;
  std::vector<std::string> invalidates;
  std::vector<const char *> read_ptrs;
  std::vector<const char *> write_ptrs;
  std::vector<const char *> invalidate_ptrs;
  std::vector<RDWR_entry> rdwr{{0, 0, 0, CCTK_VALID_GHOSTS, 0}};
  std::string where;
  std::string routine;
  std::string thorn{"TestThorn"};

  explicit FunctionFixture(const bool reverse_tags = false) {
    data.language = LangC;
    data.type = FunctionStandard;
    data.tags = make_table(reverse_tags);
    refresh();
  }

  void refresh() {
    const auto pointers = [](const std::vector<std::string> &strings,
                             std::vector<const char *> &result) {
      result.clear();
      for (const auto &value : strings)
        result.push_back(value.c_str());
    };
    pointers(reads, read_ptrs);
    pointers(writes, write_ptrs);
    pointers(invalidates, invalidate_ptrs);
    data.n_SyncGroups = static_cast<int>(sync_groups.size());
    data.SyncGroups = sync_groups.empty() ? nullptr : sync_groups.data();
    data.n_TriggerGroups = static_cast<int>(trigger_groups.size());
    data.TriggerGroups = trigger_groups.empty() ? nullptr : trigger_groups.data();
    data.n_ReadsClauses = static_cast<int>(read_ptrs.size());
    data.ReadsClauses = read_ptrs.empty() ? nullptr : read_ptrs.data();
    data.n_WritesClauses = static_cast<int>(write_ptrs.size());
    data.WritesClauses = write_ptrs.empty() ? nullptr : write_ptrs.data();
    data.n_InvalidatesClauses = static_cast<int>(invalidate_ptrs.size());
    data.InvalidatesClauses =
        invalidate_ptrs.empty() ? nullptr : invalidate_ptrs.data();
    data.n_RDWR = static_cast<int>(rdwr.size());
    data.RDWR = rdwr.empty() ? nullptr : rdwr.data();
    data.where = where.empty() ? nullptr : where.data();
    data.routine = routine.empty() ? nullptr : routine.data();
    data.thorn = thorn.empty() ? nullptr : thorn.data();
  }

  void poison() {
    where.assign("poisoned where");
    routine.assign("poisoned routine");
    thorn.assign("poisoned thorn");
    if (!reads.empty())
      reads.front().assign("poisoned read");
    if (!writes.empty())
      writes.front().assign("poisoned write");
    if (!invalidates.empty())
      invalidates.front().assign("poisoned invalidate");
    if (!sync_groups.empty())
      sync_groups.front() = 77;
    if (!trigger_groups.empty())
      trigger_groups.front() = 77;
    if (!rdwr.empty())
      rdwr.front().varindex = 77;
    auto found = tables.find(data.tags);
    if (found != tables.end() && !found->second.entries.empty())
      found->second.entries.front().bytes.assign(sizeof(int), 0xffU);
    refresh();
  }
};

struct ScopeFixture {
  std::uint64_t item_ordinal{5};
  std::uint64_t parent_local_ordinal{};
  std::string description;
  std::string implementation{"ScopeImpl"};
  FunctionFixture function;
  std::vector<int> storage_groups;
  std::vector<int> storage_timelevels;
  std::vector<int> communication_groups;
  std::vector<std::string> ifs;
  std::vector<std::string> whiles;
  std::vector<const char *> if_ptrs;
  std::vector<const char *> while_ptrs;
  cScheduleInventoryScope view{};

  explicit ScopeFixture(const std::string &root) : description(root + " scope") {
    function.where = root + "::scope_where";
    function.routine = root + "::scope_routine";
    function.refresh();
  }

  void refresh() {
    function.refresh();
    if_ptrs.clear();
    while_ptrs.clear();
    for (const auto &value : ifs)
      if_ptrs.push_back(value.c_str());
    for (const auto &value : whiles)
      while_ptrs.push_back(value.c_str());
    view = {CCTK_SCHEDULE_INVENTORY_ABI_VERSION,
            sizeof(cScheduleInventoryScope),
            item_ordinal,
            parent_local_ordinal,
            description.c_str(),
            implementation.c_str(),
            &function.data,
            static_cast<int>(storage_groups.size()),
            storage_groups.empty() ? nullptr : storage_groups.data(),
            storage_timelevels.empty() ? nullptr : storage_timelevels.data(),
            static_cast<int>(communication_groups.size()),
            communication_groups.empty() ? nullptr : communication_groups.data(),
            static_cast<int>(if_ptrs.size()),
            if_ptrs.empty() ? nullptr : if_ptrs.data(),
            static_cast<int>(while_ptrs.size()),
            while_ptrs.empty() ? nullptr : while_ptrs.data()};
  }
};

struct RecordFixture {
  std::uint64_t traversal_ordinal{};
  std::uint64_t item_ordinal{10};
  std::uint64_t parent_local_ordinal{};
  std::string root;
  std::string description;
  std::string implementation{"TestImpl"};
  FunctionFixture function;
  std::vector<int> storage_groups;
  std::vector<int> storage_timelevels;
  std::vector<int> communication_groups;
  std::vector<std::string> ifs;
  std::vector<std::string> whiles;
  std::vector<const char *> if_ptrs;
  std::vector<const char *> while_ptrs;
  std::vector<std::unique_ptr<ScopeFixture>> scopes;
  std::vector<cScheduleInventoryScope> scope_views;
  cScheduleInventoryRecord view{};
  void *opaque_handle{};

  explicit RecordFixture(std::string root_name, const bool reverse_tags = false)
      : root(std::move(root_name)), description(root + " description"),
        function(reverse_tags) {
    function.where = root + "::where";
    function.routine = root == "ODESolvers_RHS" ? "rhs_routine"
                                                  : "post_routine";
    opaque_handle = reinterpret_cast<void *>(
        static_cast<std::uintptr_t>(root == "ODESolvers_RHS" ? 0x111U : 0x222U));
    refresh();
  }

  void refresh() {
    function.refresh();
    if_ptrs.clear();
    while_ptrs.clear();
    scope_views.clear();
    for (const auto &value : ifs)
      if_ptrs.push_back(value.c_str());
    for (const auto &value : whiles)
      while_ptrs.push_back(value.c_str());
    for (auto &scope : scopes) {
      scope->refresh();
      scope_views.push_back(scope->view);
    }
    view = {CCTK_SCHEDULE_INVENTORY_ABI_VERSION,
            sizeof(cScheduleInventoryRecord),
            root.c_str(),
            0,
            traversal_ordinal,
            item_ordinal,
            parent_local_ordinal,
            description.c_str(),
            implementation.c_str(),
            &function.data,
            static_cast<int>(storage_groups.size()),
            storage_groups.empty() ? nullptr : storage_groups.data(),
            storage_timelevels.empty() ? nullptr : storage_timelevels.data(),
            static_cast<int>(communication_groups.size()),
            communication_groups.empty() ? nullptr : communication_groups.data(),
            static_cast<int>(if_ptrs.size()),
            if_ptrs.empty() ? nullptr : if_ptrs.data(),
            static_cast<int>(while_ptrs.size()),
            while_ptrs.empty() ? nullptr : while_ptrs.data(),
            scope_views.size(),
            scope_views.empty() ? nullptr : scope_views.data()};
  }

  void poison() {
    description.assign("poisoned description");
    implementation.assign("poisoned implementation");
    function.poison();
    refresh();
  }
};

struct InventoryGroup {
  int result{};
  bool poison_after_callback{};
  std::vector<std::unique_ptr<RecordFixture>> records;
};

static std::map<std::string, InventoryGroup> inventory;

static void reset() {
  if (live_clones != 0)
    throw std::runtime_error("clone leaked across test boundary");
  if (live_iterators != 0 || !iterators.empty())
    throw std::runtime_error("iterator leaked across test boundary");
  tables.clear();
  inventory.clear();
  inventory_calls.clear();
  poison_invocations = 0;
  throw_on_table_get = false;
  table_failure = TableFailure::none;
  break_group_roundtrip = false;
  break_var_roundtrip = false;
  group_infos = {{"Test::gf", {CCTK_GF, CCTK_VARIABLE_REAL, 0, 3, 1, 1, 0, 1, -1, -1}},
                 {"Test::scalar", {CCTK_SCALAR, CCTK_VARIABLE_REAL, 0, 0, 1, 1, 0, 1, -1, -1}},
                 {"Test::array", {CCTK_ARRAY, CCTK_VARIABLE_REAL, 0, 1, 1, 1, 0, 1, -1, -1}}};
  var_infos = {{"Test::u", 0, CCTK_VARIABLE_REAL, 1},
               {"Test::s", 1, CCTK_VARIABLE_REAL, 1},
               {"Test::a", 2, CCTK_VARIABLE_REAL, 1}};
}

static RecordFixture &add_body(const std::string &root,
                               const bool reverse_tags = false) {
  auto &group = inventory[root];
  group.result = 0;
  group.records.push_back(std::make_unique<RecordFixture>(root, reverse_tags));
  return *group.records.back();
}

static void setup_valid(const bool reverse_post_tags = false) {
  add_body("ODESolvers_RHS");
  auto &post = add_body("ODESolvers_PostStep", reverse_post_tags);
  post.function.data.local = 1;
  post.refresh();
}

static void destroy_source_tables() {
  for (auto it = tables.begin(); it != tables.end();) {
    if (!it->second.clone)
      it = tables.erase(it);
    else
      ++it;
  }
}

} // namespace fake

extern "C" int CCTK_ScheduleInventoryDry(const char *where,
                                           cScheduleInventoryCallback callback,
                                           void *user_data) {
  const std::string name = where == nullptr ? std::string{} : where;
  fake::inventory_calls.push_back(name);
  const auto found = fake::inventory.find(name);
  if (found == fake::inventory.end())
    return CCTK_SCHEDULE_INVENTORY_GROUP_NOT_FOUND;
  auto &group = found->second;
  if (group.result != 0)
    return group.result;
  for (auto &record : group.records) {
    record->refresh();
    const int callback_result =
        callback(&record->view, record->opaque_handle, user_data);
    if (group.poison_after_callback)
      record->poison();
    if (callback_result != 0)
      return callback_result;
  }
  return 0;
}

extern "C" int CCTK_CallFunction(void *, cFunctionData *, cGH *) {
  ++fake::poison_invocations;
  return 0;
}

extern "C" const char *CCTK_FullGroupName(const int index) {
  return index >= 0 && static_cast<std::size_t>(index) < fake::group_infos.size()
             ? fake::group_infos[static_cast<std::size_t>(index)].name.c_str()
             : nullptr;
}

extern "C" int CCTK_GroupIndex(const char *name) {
  if (name == nullptr || fake::break_group_roundtrip)
    return -1;
  for (std::size_t i = 0; i < fake::group_infos.size(); ++i)
    if (fake::group_infos[i].name == name)
      return static_cast<int>(i);
  return -1;
}

extern "C" const char *CCTK_FullVarName(const int index) {
  return index >= 0 && static_cast<std::size_t>(index) < fake::var_infos.size()
             ? fake::var_infos[static_cast<std::size_t>(index)].name.c_str()
             : nullptr;
}

extern "C" int CCTK_VarIndex(const char *name) {
  if (name == nullptr || fake::break_var_roundtrip)
    return -1;
  for (std::size_t i = 0; i < fake::var_infos.size(); ++i)
    if (fake::var_infos[i].name == name)
      return static_cast<int>(i);
  return -1;
}

extern "C" int CCTK_GroupIndexFromVarI(const int index) {
  return index >= 0 && static_cast<std::size_t>(index) < fake::var_infos.size()
             ? fake::var_infos[static_cast<std::size_t>(index)].group
             : -1;
}

extern "C" int CCTK_GroupData(const int index, cGroup *group) {
  if (group == nullptr || index < 0 ||
      static_cast<std::size_t>(index) >= fake::group_infos.size())
    return -1;
  *group = fake::group_infos[static_cast<std::size_t>(index)].data;
  return 0;
}

extern "C" int CCTK_VarTypeI(const int index) {
  return index >= 0 && static_cast<std::size_t>(index) < fake::var_infos.size()
             ? fake::var_infos[static_cast<std::size_t>(index)].type
             : -1;
}

extern "C" int CCTK_DeclaredTimeLevelsVI(const int index) {
  return index >= 0 && static_cast<std::size_t>(index) < fake::var_infos.size()
             ? fake::var_infos[static_cast<std::size_t>(index)].timelevels
             : -1;
}

extern "C" int CCTK_VarTypeSize(const int type) {
  if (type == CCTK_VARIABLE_INT)
    return static_cast<int>(sizeof(int));
  if (type == CCTK_VARIABLE_REAL)
    return static_cast<int>(sizeof(double));
  if (type == CCTK_VARIABLE_POINTER || type == CCTK_VARIABLE_POINTER_TO_CONST ||
      type == CCTK_VARIABLE_FPOINTER)
    return static_cast<int>(sizeof(void *));
  return -1;
}

extern "C" int Util_TableQueryFlags(const int handle) {
  const auto found = fake::tables.find(handle);
  return found == fake::tables.end() ? -1 : found->second.flags;
}

extern "C" int Util_TableQueryNKeys(const int handle) {
  const auto found = fake::tables.find(handle);
  return found == fake::tables.end()
             ? -1
             : static_cast<int>(found->second.entries.size());
}

extern "C" int Util_TableQueryMaxKeyLength(const int handle) {
  const auto found = fake::tables.find(handle);
  if (found == fake::tables.end())
    return -1;
  std::size_t length = 0;
  for (const auto &entry : found->second.entries)
    length = std::max(length, entry.key.size());
  return static_cast<int>(length);
}

extern "C" int Util_TableItCreate(const int handle) {
  if (fake::table_failure == fake::TableFailure::iterator_create)
    return -1;
  if (fake::tables.find(handle) == fake::tables.end())
    return -1;
  const int iterator = fake::next_iterator_handle++;
  fake::iterators.emplace(iterator, fake::Iterator{handle, 0});
  ++fake::live_iterators;
  return iterator;
}

extern "C" int Util_TableItDestroy(const int handle) {
  if (fake::iterators.erase(handle) != 1)
    return -1;
  --fake::live_iterators;
  return 0;
}

extern "C" int Util_TableItQueryKeyValueInfo(const int iterator_handle,
                                               const int capacity, char *key,
                                               CCTK_INT *type,
                                               CCTK_INT *count) {
  if (fake::table_failure == fake::TableFailure::iterator_query)
    return -1;
  const auto iterator = fake::iterators.find(iterator_handle);
  if (iterator == fake::iterators.end())
    return -1;
  const auto table = fake::tables.find(iterator->second.table);
  if (table == fake::tables.end() ||
      iterator->second.index >= table->second.entries.size())
    return -1;
  const auto &entry = table->second.entries[iterator->second.index];
  if (capacity <= static_cast<int>(entry.key.size()) || key == nullptr ||
      type == nullptr || count == nullptr)
    return -1;
  std::memcpy(key, entry.key.c_str(), entry.key.size() + 1U);
  *type = entry.type;
  *count = entry.count;
  return static_cast<int>(entry.key.size());
}

extern "C" int Util_TableGetGenericArray(const int handle, const int type,
                                          const int count, void *output,
                                          const char *key) {
  if (fake::throw_on_table_get) {
    fake::throw_on_table_get = false;
    throw std::runtime_error("injected table read exception");
  }
  const auto table = fake::tables.find(handle);
  if (table == fake::tables.end() || key == nullptr)
    return -1;
  const auto lower = [](std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](const char c) {
      return c >= 'A' && c <= 'Z' ? static_cast<char>(c - 'A' + 'a') : c;
    });
    return value;
  };
  const std::string wanted = lower(key);
  for (const auto &entry : table->second.entries) {
    if (lower(entry.key) != wanted)
      continue;
    if (entry.get_result.has_value())
      return *entry.get_result;
    if (entry.type != type || entry.count != count || output == nullptr)
      return -1;
    std::memcpy(output, entry.bytes.data(), entry.bytes.size());
    return entry.count;
  }
  return -1;
}

extern "C" int Util_TableItAdvance(const int handle) {
  if (fake::table_failure == fake::TableFailure::iterator_advance)
    return -1;
  const auto iterator = fake::iterators.find(handle);
  if (iterator == fake::iterators.end())
    return -1;
  const auto table = fake::tables.find(iterator->second.table);
  if (table == fake::tables.end())
    return -1;
  ++iterator->second.index;
  return iterator->second.index < table->second.entries.size() ? 1 : 0;
}

extern "C" int Util_TableItQueryIsNull(const int handle) {
  if (fake::table_failure == fake::TableFailure::iterator_end_state)
    return 0;
  const auto iterator = fake::iterators.find(handle);
  if (iterator == fake::iterators.end())
    return -1;
  const auto table = fake::tables.find(iterator->second.table);
  return table != fake::tables.end() &&
                 iterator->second.index >= table->second.entries.size()
             ? 1
             : 0;
}

extern "C" int Util_TableClone(const int handle) {
  if (fake::table_failure == fake::TableFailure::clone)
    return -1;
  const auto found = fake::tables.find(handle);
  if (found == fake::tables.end())
    return -1;
  const int clone = fake::next_table_handle++;
  auto value = found->second;
  value.clone = true;
  fake::tables.emplace(clone, std::move(value));
  ++fake::live_clones;
  return clone;
}

extern "C" int Util_TableDestroy(const int handle) {
  const auto found = fake::tables.find(handle);
  if (found == fake::tables.end())
    return -1;
  if (found->second.clone)
    --fake::live_clones;
  fake::tables.erase(found);
  return 0;
}

namespace {

class TestFailure : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

void check(const bool condition, const char *expression, const int line) {
  if (!condition)
    throw TestFailure(std::string("line ") + std::to_string(line) +
                      ": " + expression);
}

#define CHECK(expression)                                                       \
  check(static_cast<bool>(expression), #expression, __LINE__)

CarpetX::ScheduleBuildProvenance valid_provenance() {
  return {1,
          CCTK_SCHEDULE_INVENTORY_ABI_VERSION,
          std::string(40, 'a'),
          std::string(64, 'b'),
          std::string(40, 'c'),
          {{"8B3A", std::string(64, 'd')}},
          "GCC",
          "13.2.0",
          "x86_64-linux-gnu",
          "itanium",
          "little",
          {{CCTK_VARIABLE_INT, static_cast<int>(sizeof(int))},
           {CCTK_VARIABLE_REAL, static_cast<int>(sizeof(double))}},
          std::string(64, 'e')};
}

std::vector<CarpetX::CanonicalTag> expected_tags() {
  const double real_value = 1.5;
  const int integer_value = 9;
  return {{"alpha", CCTK_VARIABLE_REAL, 1, fake::bytes_of(real_value)},
          {"beta", CCTK_VARIABLE_INT, 1, fake::bytes_of(integer_value)}};
}

CarpetX::CanonicalScheduleItem expected_item(const std::string &root,
                                             const int local,
                                             const bool scope = false) {
  CarpetX::CanonicalScheduleItem item;
  item.description = root + (scope ? " scope" : " description");
  item.implementation = scope ? "ScopeImpl" : "TestImpl";
  item.where = root + (scope ? "::scope_where" : "::where");
  item.thorn = "TestThorn";
  item.routine = scope ? root + "::scope_routine"
                       : (root == "ODESolvers_RHS" ? "rhs_routine"
                                                    : "post_routine");
  item.language = LangC;
  item.function_type = FunctionStandard;
  item.execution_flags.local = local;
  item.effective_mode = local == 0 ? CarpetX::CanonicalScheduleMode::default_local
                                   : CarpetX::CanonicalScheduleMode::local;
  item.tags = expected_tags();
  item.reads_clauses = {"Test::u"};
  item.accesses = {{"Test::u", 0, CCTK_VALID_GHOSTS, 0, 0}};
  return item;
}

CarpetX::ScheduleCertificationExpectation valid_expectation() {
  CarpetX::ScheduleCertificationExpectation expectation;
  expectation.provenance = valid_provenance();
  constexpr const char *roots[] = {"ODESolvers_RHS", "ODESolvers_PostStep"};
  for (std::size_t i = 0; i < 2; ++i) {
    auto &target = expectation.manifest.targets[i];
    target.target = i == 0 ? CarpetX::SubcyclingScheduleTarget::rhs
                           : CarpetX::SubcyclingScheduleTarget::post_step;
    target.body_root = roots[i];
    target.entry_state = CarpetX::SpecialScheduleGroupState::missing;
    target.exit_state = CarpetX::SpecialScheduleGroupState::missing;
    target.entry_exit_included = false;
    target.records.push_back({0, 10, 0,
                              expected_item(roots[i], i == 0 ? 0 : 1), {}});
  }
  return expectation;
}

CarpetX::ScheduleCertificationResult certify(
    const CarpetX::ScheduleCertificationExpectation &expectation) {
  return CarpetX::certify_local_subcycling_schedules(expectation,
                                                      valid_provenance());
}

void check_failure(const CarpetX::ScheduleCertificationResult &result,
                   const CarpetX::ScheduleCertificationErrorCode code) {
  CHECK(!result);
  CHECK(result.registry == nullptr);
  CHECK(result.failure.has_value());
  CHECK(result.failure->code == code);
}

void test_exact_six_call_order_and_two_body_targets() {
  fake::setup_valid();
  const auto result = certify(valid_expectation());
  CHECK(result);
  CHECK(result.registry->size(CarpetX::SubcyclingScheduleTarget::rhs) == 1);
  CHECK(result.registry->size(CarpetX::SubcyclingScheduleTarget::post_step) == 1);
  const std::vector<std::string> expected{
      "ODESolvers_RHS$ENTRY", "ODESolvers_RHS", "ODESolvers_RHS$EXIT",
      "ODESolvers_PostStep$ENTRY", "ODESolvers_PostStep",
      "ODESolvers_PostStep$EXIT"};
  CHECK(fake::inventory_calls == expected);
}

void test_missing_entry_and_exit_are_required() {
  fake::setup_valid();
  const auto result = certify(valid_expectation());
  CHECK(result);
  for (const auto &target : result.registry->manifest().targets) {
    CHECK(target.entry_state == CarpetX::SpecialScheduleGroupState::missing);
    CHECK(target.exit_state == CarpetX::SpecialScheduleGroupState::missing);
    CHECK(!target.entry_exit_included);
  }
}

void test_present_empty_entry_or_exit_rejects() {
  {
    fake::setup_valid();
    fake::inventory["ODESolvers_RHS$ENTRY"].result = 0;
    const auto result = certify(valid_expectation());
    check_failure(result,
                  CarpetX::ScheduleCertificationErrorCode::special_group_present);
    CHECK(fake::live_clones == 0);
  }
  fake::reset();
  {
    fake::setup_valid();
    fake::inventory["ODESolvers_RHS$EXIT"].result = 0;
    const auto result = certify(valid_expectation());
    check_failure(result,
                  CarpetX::ScheduleCertificationErrorCode::special_group_present);
    CHECK(fake::live_clones == 0);
  }
}

void test_body_missing_or_error_rejects() {
  {
    fake::setup_valid();
    fake::inventory.erase("ODESolvers_RHS");
    const auto result = certify(valid_expectation());
    check_failure(result, CarpetX::ScheduleCertificationErrorCode::inventory_error);
  }
  fake::reset();
  {
    fake::setup_valid();
    fake::inventory["ODESolvers_RHS"].result = -7;
    const auto result = certify(valid_expectation());
    check_failure(result, CarpetX::ScheduleCertificationErrorCode::inventory_error);
  }
}

void test_borrowed_metadata_is_deep_copied() {
  fake::setup_valid();
  fake::inventory["ODESolvers_RHS"].poison_after_callback = true;
  fake::inventory["ODESolvers_PostStep"].poison_after_callback = true;
  const auto result = certify(valid_expectation());
  CHECK(result);
  CHECK(result.registry->manifest() == valid_expectation().manifest);
}

void test_owned_function_data_rebinds_every_pointer() {
  fake::RecordFixture record("ODESolvers_RHS");
  record.function.sync_groups = {0};
  record.function.trigger_groups = {0};
  record.function.writes = {"write clause"};
  record.function.invalidates = {"invalidate clause"};
  record.storage_groups = {0};
  record.storage_timelevels = {0};
  record.communication_groups = {0};
  record.ifs = {"if clause"};
  record.whiles = {"while clause"};
  record.refresh();
  const auto *const source_sync = record.function.data.SyncGroups;
  const auto *const source_trigger = record.function.data.TriggerGroups;
  const auto *const source_rdwr = record.function.data.RDWR;
  const auto *const source_reads = record.function.data.ReadsClauses;
  const auto *const source_writes = record.function.data.WritesClauses;
  const auto *const source_invalidates = record.function.data.InvalidatesClauses;
  const auto *const source_where = record.function.data.where;
  const auto *const source_routine = record.function.data.routine;
  const auto *const source_thorn = record.function.data.thorn;
  const int clone = Util_TableClone(record.function.data.tags);
  CHECK(clone >= 0);
  {
    CarpetX::OwnedFunctionData owned;
    owned.capture(record.function.data, CarpetX::TableHandle(clone),
                  record.opaque_handle, record.view);
    CHECK(owned.data.SyncGroups != source_sync);
    CHECK(owned.data.TriggerGroups != source_trigger);
    CHECK(owned.data.RDWR != source_rdwr);
    CHECK(owned.data.ReadsClauses != source_reads);
    CHECK(owned.data.WritesClauses != source_writes);
    CHECK(owned.data.InvalidatesClauses != source_invalidates);
    CHECK(owned.data.ReadsClauses[0] != source_reads[0]);
    CHECK(owned.data.WritesClauses[0] != source_writes[0]);
    CHECK(owned.data.InvalidatesClauses[0] != source_invalidates[0]);
    CHECK(owned.data.where != source_where);
    CHECK(owned.data.routine != source_routine);
    CHECK(owned.data.thorn != source_thorn);
    CHECK(owned.data.tags == clone);
    CHECK(owned.opaque_call_handle == record.opaque_handle);
    CHECK(owned.storage_group_indices == std::vector<int>{0});
    CHECK(owned.storage_timelevels == std::vector<int>{0});
    CHECK(owned.communication_group_indices == std::vector<int>{0});
    CHECK(owned.if_strings == std::vector<std::string>{"if clause"});
    CHECK(owned.while_strings == std::vector<std::string>{"while clause"});
    record.function.poison();
    CHECK(std::string(owned.data.where) == "ODESolvers_RHS::where");
    CHECK(std::string(owned.data.routine) == "rhs_routine");
    CHECK(std::string(owned.data.thorn) == "TestThorn");
    CHECK(std::string(owned.data.ReadsClauses[0]) == "Test::u");
    CHECK(std::string(owned.data.WritesClauses[0]) == "write clause");
    CHECK(std::string(owned.data.InvalidatesClauses[0]) == "invalidate clause");
    CHECK(owned.data.RDWR[0].varindex == 0);
  }
  CHECK(fake::live_clones == 0);
}

void test_leaf_tag_clone_and_scope_tag_copy_lifetimes() {
  fake::setup_valid();
  for (const char *root : {"ODESolvers_RHS", "ODESolvers_PostStep"}) {
    auto &record = *fake::inventory[root].records.front();
    record.scopes.push_back(std::make_unique<fake::ScopeFixture>(root));
    record.refresh();
  }
  auto expectation = valid_expectation();
  for (std::size_t i = 0; i < 2; ++i) {
    const std::string root = i == 0 ? "ODESolvers_RHS" : "ODESolvers_PostStep";
    expectation.manifest.targets[i].records[0].ancestry.push_back(
        {5, 0, expected_item(root, 0, true)});
  }
  {
    const auto result = certify(expectation);
    CHECK(result);
    CHECK(fake::live_clones == 2);
    fake::destroy_source_tables();
    CHECK(result.registry->manifest() == expectation.manifest);
    CHECK(fake::live_clones == 2);
  }
  CHECK(fake::live_clones == 0);
}

CarpetX::ScheduleCertificationExpectation setup_two_record_rhs_ancestry() {
  fake::setup_valid();
  auto &first = *fake::inventory["ODESolvers_RHS"].records.front();
  first.scopes.push_back(
      std::make_unique<fake::ScopeFixture>("ODESolvers_RHS"));
  first.refresh();

  auto &second = fake::add_body("ODESolvers_RHS");
  second.traversal_ordinal = 1;
  second.item_ordinal = 11;
  second.parent_local_ordinal = 1;
  second.scopes.push_back(
      std::make_unique<fake::ScopeFixture>("ODESolvers_RHS"));
  second.refresh();

  auto expectation = valid_expectation();
  expectation.manifest.targets[0].records[0].ancestry.push_back(
      {5, 0, expected_item("ODESolvers_RHS", 0, true)});
  expectation.manifest.targets[0].records.push_back(
      {1,
       11,
       1,
       expected_item("ODESolvers_RHS", 0),
       {{5, 0, expected_item("ODESolvers_RHS", 0, true)}}});
  return expectation;
}

void test_multi_record_ancestry_uses_sibling_local_ordinals() {
  const auto result = certify(setup_two_record_rhs_ancestry());
  CHECK(result);
  CHECK(result.registry->manifest().targets[0].records[1]
            .parent_local_ordinal == 1);
  CHECK(result.registry->manifest().targets[0].records[1]
            .ancestry[0]
            .parent_local_ordinal == 0);
}

void test_contradictory_multi_record_ancestry_rejects_exact_fields() {
  const std::vector<std::string> expected_fields{
      "record.ancestry[0].item_identity",
      "record.ancestry[0].parent_path",
      "record.ancestry[0].parent_local_ordinal",
      "record.ancestry[0].parent_local_ordinal",
      "record.ancestry[0].parent_local_ordinal"};
  for (int which = 0; which < 5; ++which) {
    if (which != 0)
      fake::reset();
    auto expectation = setup_two_record_rhs_ancestry();
    auto &second = *fake::inventory["ODESolvers_RHS"].records[1];
    auto &expected_second = expectation.manifest.targets[0].records[1];
    if (which == 0) {
      second.scopes[0]->function.routine = "contradictory_scope_routine";
      expected_second.ancestry[0].item.routine =
          "contradictory_scope_routine";
    } else if (which == 1) {
      auto &first = *fake::inventory["ODESolvers_RHS"].records[0];
      first.scopes.clear();
      auto first_outer =
          std::make_unique<fake::ScopeFixture>("ODESolvers_RHS");
      first_outer->item_ordinal = 4;
      first.scopes.push_back(std::move(first_outer));
      first.scopes.push_back(
          std::make_unique<fake::ScopeFixture>("ODESolvers_RHS"));
      first.refresh();
      second.parent_local_ordinal = 0;
      second.scopes.clear();
      auto repeated =
          std::make_unique<fake::ScopeFixture>("ODESolvers_RHS");
      repeated->item_ordinal = 5;
      repeated->parent_local_ordinal = 1;
      second.scopes.push_back(std::move(repeated));
      expectation.manifest.targets[0].records[0].ancestry = {
          {4, 0, expected_item("ODESolvers_RHS", 0, true)},
          {5, 0, expected_item("ODESolvers_RHS", 0, true)}};
      expected_second.parent_local_ordinal = 0;
      expected_second.ancestry = {
          {5, 1, expected_item("ODESolvers_RHS", 0, true)}};
    } else if (which == 2) {
      second.scopes[0]->parent_local_ordinal = 1;
      expected_second.ancestry[0].parent_local_ordinal = 1;
    } else {
      second.item_ordinal = 21;
      second.scopes[0]->item_ordinal = 11;
      second.scopes[0]->parent_local_ordinal = which == 3 ? 0 : 2;
      expected_second.item_ordinal = 21;
      expected_second.ancestry[0].item_ordinal = 11;
      expected_second.ancestry[0].parent_local_ordinal = which == 3 ? 0 : 2;
    }
    second.refresh();
    const auto result = certify(expectation);
    check_failure(result,
                  CarpetX::ScheduleCertificationErrorCode::malformed_record);
    CHECK(result.failure->target == CarpetX::SubcyclingScheduleTarget::rhs);
    CHECK(result.failure->traversal_ordinal == 1);
    CHECK(result.failure->field ==
          expected_fields[static_cast<std::size_t>(which)]);
  }
}

void test_repeated_item_cannot_change_from_leaf_to_scope() {
  fake::setup_valid();
  auto &second = fake::add_body("ODESolvers_RHS");
  second.traversal_ordinal = 1;
  second.item_ordinal = 11;
  second.parent_local_ordinal = 0;
  auto replayed_leaf =
      std::make_unique<fake::ScopeFixture>("ODESolvers_RHS");
  replayed_leaf->item_ordinal = 10;
  replayed_leaf->parent_local_ordinal = 0;
  replayed_leaf->description = "ODESolvers_RHS description";
  replayed_leaf->implementation = "TestImpl";
  replayed_leaf->function.where = "ODESolvers_RHS::where";
  replayed_leaf->function.routine = "rhs_routine";
  second.scopes.push_back(std::move(replayed_leaf));
  second.refresh();

  auto expectation = valid_expectation();
  expectation.manifest.targets[0].records.push_back(
      {1,
       11,
       0,
       expected_item("ODESolvers_RHS", 0),
       {{10, 0, expected_item("ODESolvers_RHS", 0)}}});
  const auto result = certify(expectation);
  check_failure(result,
                CarpetX::ScheduleCertificationErrorCode::malformed_record);
  CHECK(result.failure->target == CarpetX::SubcyclingScheduleTarget::rhs);
  CHECK(result.failure->traversal_ordinal == 1);
  CHECK(result.failure->field == "record.ancestry[0].node_kind");
}

void test_new_scope_item_ordinal_cannot_move_backwards_across_records() {
  fake::setup_valid();
  auto &second = fake::add_body("ODESolvers_RHS");
  second.traversal_ordinal = 1;
  second.item_ordinal = 11;
  second.parent_local_ordinal = 0;
  second.scopes.push_back(
      std::make_unique<fake::ScopeFixture>("ODESolvers_RHS"));
  second.scopes[0]->parent_local_ordinal = 1;
  second.refresh();

  auto expectation = valid_expectation();
  expectation.manifest.targets[0].records.push_back(
      {1,
       11,
       0,
       expected_item("ODESolvers_RHS", 0),
       {{5, 1, expected_item("ODESolvers_RHS", 0, true)}}});
  const auto result = certify(expectation);
  check_failure(result,
                CarpetX::ScheduleCertificationErrorCode::malformed_record);
  CHECK(result.failure->target == CarpetX::SubcyclingScheduleTarget::rhs);
  CHECK(result.failure->traversal_ordinal == 1);
  CHECK(result.failure->field == "record.ancestry[0].item_ordinal");
}

void test_table_boundary_failures_release_iterators_and_clones() {
  const std::vector<std::pair<fake::TableFailure, std::string>> cases{
      {fake::TableFailure::iterator_create, "record.item.tags.iterator"},
      {fake::TableFailure::iterator_query, "record.item.tags.iterator"},
      {fake::TableFailure::iterator_advance, "record.item.tags.iterator"},
      {fake::TableFailure::iterator_end_state, "record.item.tags.iterator"},
      {fake::TableFailure::clone, "record.item.tags.clone"}};
  for (std::size_t i = 0; i < cases.size(); ++i) {
    if (i != 0)
      fake::reset();
    fake::setup_valid();
    fake::table_failure = cases[i].first;
    const auto result = certify(valid_expectation());
    check_failure(result, CarpetX::ScheduleCertificationErrorCode::tag_error);
    CHECK(result.failure->field == cases[i].second);
    CHECK(fake::live_iterators == 0);
    CHECK(fake::iterators.empty());
    CHECK(fake::live_clones == 0);
  }
}

void test_late_second_target_failure_discards_first_target() {
  fake::setup_valid();
  auto &post = *fake::inventory["ODESolvers_PostStep"].records.front();
  fake::tables.at(post.function.data.tags).flags = 0;
  const auto result = certify(valid_expectation());
  check_failure(result, CarpetX::ScheduleCertificationErrorCode::tag_error);
  CHECK(fake::live_clones == 0);
  CHECK(fake::inventory_calls.size() == 5);
}

void test_poison_handles_are_never_invoked() {
  fake::setup_valid();
  fake::inventory["ODESolvers_RHS"].records.front()->opaque_handle =
      reinterpret_cast<void *>(static_cast<std::uintptr_t>(0xdeadU));
  fake::inventory["ODESolvers_PostStep"].records.front()->opaque_handle =
      reinterpret_cast<void *>(static_cast<std::uintptr_t>(0xbeefU));
  const auto result = certify(valid_expectation());
  CHECK(result);
  CHECK(fake::poison_invocations == 0);
}

void test_default_local_and_explicit_local_c_gf_real_tl0_pass() {
  fake::setup_valid();
  const auto result = certify(valid_expectation());
  CHECK(result);
  const auto &targets = result.registry->manifest().targets;
  CHECK(targets[0].records[0].item.effective_mode ==
        CarpetX::CanonicalScheduleMode::default_local);
  CHECK(targets[1].records[0].item.effective_mode ==
        CarpetX::CanonicalScheduleMode::local);
  CHECK(targets[0].records[0].item.accesses[0].variable_name == "Test::u");
}

void test_random_tag_iteration_order_has_equal_manifest() {
  fake::setup_valid(true);
  const auto result = certify(valid_expectation());
  CHECK(result);
  CHECK(result.registry->manifest().targets[0].records[0].item.tags ==
        result.registry->manifest().targets[1].records[0].item.tags);
}

void test_tag_flags_types_counts_values_and_duplicates_fail_closed() {
  const auto run_case = [](const int which) {
    fake::setup_valid();
    auto &record = *fake::inventory["ODESolvers_RHS"].records.front();
    auto &table = fake::tables.at(record.function.data.tags);
    if (which == 0)
      table.flags = 0;
    else if (which == 1)
      table.entries[0].type = CCTK_VARIABLE_POINTER;
    else if (which == 2)
      table.entries[0].count = 0;
    else if (which == 3)
      table.entries[0].get_result = 0;
    else
      table.entries = {{"Dup", CCTK_VARIABLE_INT, 1,
                        fake::bytes_of(1), {}},
                       {"dUP", CCTK_VARIABLE_INT, 1,
                        fake::bytes_of(2), {}}};
    const auto result = certify(valid_expectation());
    check_failure(result, CarpetX::ScheduleCertificationErrorCode::tag_error);
    CHECK(fake::live_clones == 0);
  };
  for (int which = 0; which < 5; ++which) {
    if (which != 0)
      fake::reset();
    run_case(which);
  }
}

void test_group_and_variable_full_name_roundtrip_is_required() {
  {
    fake::setup_valid();
    fake::break_var_roundtrip = true;
    const auto result = certify(valid_expectation());
    check_failure(
        result, CarpetX::ScheduleCertificationErrorCode::name_resolution_error);
  }
  fake::reset();
  {
    fake::setup_valid();
    fake::break_group_roundtrip = true;
    const auto result = certify(valid_expectation());
    check_failure(
        result, CarpetX::ScheduleCertificationErrorCode::name_resolution_error);
  }
}

void test_scalar_array_tl1_and_unknown_region_bits_reject() {
  for (int which = 0; which < 4; ++which) {
    if (which != 0)
      fake::reset();
    fake::setup_valid();
    auto &record = *fake::inventory["ODESolvers_RHS"].records.front();
    if (which == 0)
      record.function.rdwr[0].varindex = 1;
    else if (which == 1)
      record.function.rdwr[0].varindex = 2;
    else if (which == 2)
      record.function.rdwr[0].timelevel = 1;
    else
      record.function.rdwr[0].where_rd = 8;
    record.refresh();
    const auto result = certify(valid_expectation());
    check_failure(
        result, CarpetX::ScheduleCertificationErrorCode::unsupported_metadata);
  }
}

void test_nonlocal_early_late_singlemap_and_loop_modes_reject() {
  using FlagCase = std::pair<int cFunctionData::*, int>;
  const std::vector<FlagCase> flags{{&cFunctionData::local, 2},
                                    {&cFunctionData::meta, 1},
                                    {&cFunctionData::meta_early, 1},
                                    {&cFunctionData::meta_late, 1},
                                    {&cFunctionData::global, 1},
                                    {&cFunctionData::global_early, 1},
                                    {&cFunctionData::global_late, 1},
                                    {&cFunctionData::level, 1},
                                    {&cFunctionData::singlemap, 1},
                                    {&cFunctionData::loop_meta, 1},
                                    {&cFunctionData::loop_global, 1},
                                    {&cFunctionData::loop_level, 1},
                                    {&cFunctionData::loop_singlemap, 1},
                                    {&cFunctionData::loop_local, 1}};
  for (std::size_t i = 0; i < flags.size(); ++i) {
    if (i != 0)
      fake::reset();
    fake::setup_valid();
    auto &data = fake::inventory["ODESolvers_RHS"].records.front()->function.data;
    data.*flags[i].first = flags[i].second;
    const auto result = certify(valid_expectation());
    check_failure(
        result, CarpetX::ScheduleCertificationErrorCode::unsupported_metadata);
  }
}

void test_if_while_storage_communication_sync_and_trigger_reject() {
  for (int which = 0; which < 6; ++which) {
    if (which != 0)
      fake::reset();
    fake::setup_valid();
    auto &record = *fake::inventory["ODESolvers_RHS"].records.front();
    if (which == 0)
      record.ifs = {"condition"};
    else if (which == 1)
      record.whiles = {"condition"};
    else if (which == 2) {
      record.storage_groups = {0};
      record.storage_timelevels = {0};
    } else if (which == 3)
      record.communication_groups = {0};
    else if (which == 4)
      record.function.sync_groups = {0};
    else
      record.function.trigger_groups = {0};
    record.refresh();
    const auto result = certify(valid_expectation());
    check_failure(
        result, CarpetX::ScheduleCertificationErrorCode::unsupported_metadata);
  }
}

void test_wavetoy_global_sync_is_explicitly_rejected() {
  fake::setup_valid();
  auto &function =
      fake::inventory["ODESolvers_RHS"].records.front()->function;
  function.thorn = "WaveToyX";
  function.routine = "WaveToyX_Boundaries";
  function.sync_groups = {0};
  function.data.global = 1;
  function.refresh();

  const auto result = certify(valid_expectation());
  check_failure(
      result, CarpetX::ScheduleCertificationErrorCode::unsupported_metadata);
  CHECK(result.failure->target ==
        CarpetX::SubcyclingScheduleTarget::rhs);
  CHECK(result.failure->traversal_ordinal == 0);
  CHECK(result.failure->field == "record.item.execution_mode");
  CHECK(result.failure->detail ==
        "only default-local or explicit local is supported");
  CHECK(result.registry == nullptr);
  CHECK(fake::poison_invocation_count == 0);
  CHECK(fake::live_cloned_table_count == 0);
}

void test_raw_and_normalized_access_contradiction_rejects() {
  for (int which = 0; which < 2; ++which) {
    if (which != 0)
      fake::reset();
    fake::setup_valid();
    auto &function =
        fake::inventory["ODESolvers_RHS"].records.front()->function;
    if (which == 0)
      function.reads.clear();
    else
      function.rdwr[0].where_rd = 0;
    function.refresh();
    const auto result = certify(valid_expectation());
    check_failure(result,
                  CarpetX::ScheduleCertificationErrorCode::malformed_record);
  }
}

void test_extra_missing_reordered_or_changed_manifest_rejects() {
  const std::vector<std::string> expected_fields{
      "manifest.targets[0].records.size",
      "manifest.targets[0].records.size",
      "manifest.targets[0].target",
      "manifest.targets[1].records[0].item.routine"};
  for (int which = 0; which < 4; ++which) {
    if (which != 0)
      fake::reset();
    fake::setup_valid();
    auto expectation = valid_expectation();
    if (which == 0)
      expectation.manifest.targets[0].records.push_back(
          expectation.manifest.targets[0].records[0]);
    else if (which == 1)
      expectation.manifest.targets[0].records.clear();
    else if (which == 2)
      std::swap(expectation.manifest.targets[0],
                expectation.manifest.targets[1]);
    else
      expectation.manifest.targets[1].records[0].item.routine = "changed";
    const auto result = certify(expectation);
    check_failure(result,
                  CarpetX::ScheduleCertificationErrorCode::manifest_mismatch);
    CHECK(result.failure->field ==
          expected_fields[static_cast<std::size_t>(which)]);
    CHECK(fake::live_clones == 0);
  }
}

void test_pointer_and_table_handle_values_do_not_affect_equality() {
  fake::setup_valid();
  const auto expectation = valid_expectation();
  CarpetX::CanonicalScheduleBundle first_manifest;
  {
    const auto first = certify(expectation);
    CHECK(first);
    first_manifest = first.registry->manifest();
  }
  CHECK(fake::live_clones == 0);
  fake::reset();
  fake::setup_valid(true);
  fake::inventory["ODESolvers_RHS"].records.front()->opaque_handle =
      reinterpret_cast<void *>(static_cast<std::uintptr_t>(0x123456U));
  fake::inventory["ODESolvers_PostStep"].records.front()->opaque_handle =
      reinterpret_cast<void *>(static_cast<std::uintptr_t>(0x654321U));
  const auto second = certify(expectation);
  CHECK(second);
  CHECK(first_manifest == second.registry->manifest());
}

void test_invalid_or_mismatched_build_provenance_rejects() {
  {
    fake::setup_valid();
    auto expectation = valid_expectation();
    expectation.provenance.executable_sha256 = "bad";
    const auto result = CarpetX::certify_local_subcycling_schedules(
        expectation, valid_provenance());
    check_failure(
        result, CarpetX::ScheduleCertificationErrorCode::provenance_mismatch);
    CHECK(result.failure->field == "provenance.executable_sha256");
    CHECK(fake::inventory_calls.empty());
  }
  fake::reset();
  {
    fake::setup_valid();
    auto observed = valid_provenance();
    observed.compiler_version = "different";
    const auto result = CarpetX::certify_local_subcycling_schedules(
        valid_expectation(), observed);
    check_failure(
        result, CarpetX::ScheduleCertificationErrorCode::provenance_mismatch);
    CHECK(result.failure->field == "provenance.compiler_version");
    CHECK(fake::inventory_calls.empty());
  }
  fake::reset();
  {
    fake::setup_valid();
    auto expectation = valid_expectation();
    std::transform(expectation.provenance.cactus_flesh_parent.begin(),
                   expectation.provenance.cactus_flesh_parent.end(),
                   expectation.provenance.cactus_flesh_parent.begin(),
                   [](const char) { return 'A'; });
    const auto result = CarpetX::certify_local_subcycling_schedules(
        expectation, valid_provenance());
    CHECK(result);
  }
}

void test_provenance_type_sizes_match_live_build_and_cover_observed_tags() {
  {
    fake::setup_valid();
    auto expectation = valid_expectation();
    expectation.provenance.cctk_type_sizes[0].second += 1;
    const auto result = CarpetX::certify_local_subcycling_schedules(
        expectation, expectation.provenance);
    check_failure(
        result, CarpetX::ScheduleCertificationErrorCode::provenance_mismatch);
    CHECK(result.failure->field == "provenance.cctk_type_sizes[0].size");
    CHECK(fake::inventory_calls.empty());
  }
  fake::reset();
  {
    fake::setup_valid();
    auto expectation = valid_expectation();
    expectation.provenance.cctk_type_sizes.pop_back();
    const auto result = CarpetX::certify_local_subcycling_schedules(
        expectation, expectation.provenance);
    check_failure(
        result, CarpetX::ScheduleCertificationErrorCode::provenance_mismatch);
    CHECK(result.failure->field ==
          "provenance.cctk_type_sizes.missing[130]");
    CHECK(fake::inventory_calls.size() == 6);
  }
}

void test_registry_publishes_only_after_both_targets_match() {
  {
    fake::setup_valid();
    auto expectation = valid_expectation();
    expectation.manifest.targets[1].records[0].item.routine = "wrong";
    const auto result = certify(expectation);
    check_failure(result,
                  CarpetX::ScheduleCertificationErrorCode::manifest_mismatch);
    CHECK(fake::live_clones == 0);
  }
  fake::reset();
  {
    fake::setup_valid();
    const auto result = certify(valid_expectation());
    CHECK(result);
    CHECK(result.registry != nullptr);
    CHECK(result.registry->size(CarpetX::SubcyclingScheduleTarget::rhs) == 1);
    CHECK(result.registry->size(CarpetX::SubcyclingScheduleTarget::post_step) == 1);
  }
}

void test_native_gate_observation_returns_only_canonical_expectation() {
  fake::setup_valid();
  const auto observed =
      CarpetX::observe_local_subcycling_schedules_for_native_gate(
          valid_provenance());
  CHECK(observed);
  CHECK(!observed.failure.has_value());
  CHECK(observed.expectation.has_value());
  CHECK(observed.expectation.value().manifest == valid_expectation().manifest);
  CHECK(observed.expectation.value().provenance == valid_provenance());
  CHECK(fake::inventory_calls.size() == 6);
  CHECK(fake::live_clones == 0);
}

void test_callback_exception_is_contained_and_rolls_back() {
  fake::setup_valid();
  fake::throw_on_table_get = true;
  const auto result = certify(valid_expectation());
  check_failure(
      result,
      CarpetX::ScheduleCertificationErrorCode::callback_internal_error);
  CHECK(result.failure->field == "record.exception");
  CHECK(fake::live_clones == 0);
}

using TestFunction = void (*)();

int run_test(const char *name, TestFunction test) {
  std::string failure_message;
  try {
    fake::reset();
    test();
  } catch (const std::exception &error) {
    failure_message = error.what();
  }
  const auto record_failure = [&](const std::string &message) {
    if (!failure_message.empty())
      failure_message += "; ";
    failure_message += message;
  };
  if (fake::poison_invocation_count != 0)
    record_failure("opaque schedule handle was invoked");
  if (fake::live_cloned_table_count != 0)
    record_failure("owned tag clone leaked past test lifetime");
  if (fake::live_iterators != 0 || !fake::iterators.empty())
    record_failure("tag-table iterator leaked past test lifetime");
  if (!failure_message.empty()) {
    std::cerr << "[FAIL] " << name << ": " << failure_message << '\n';
    return 1;
  }
  std::cout << "[PASS] " << name << '\n';
  return 0;
}

} // namespace

int main() {
  const std::vector<std::pair<const char *, TestFunction>> tests{
      {"exact six-call order and two body targets",
       test_exact_six_call_order_and_two_body_targets},
      {"missing entry and exit are required",
       test_missing_entry_and_exit_are_required},
      {"present empty entry or exit rejects",
       test_present_empty_entry_or_exit_rejects},
      {"body missing or error rejects", test_body_missing_or_error_rejects},
      {"borrowed metadata is deep copied",
       test_borrowed_metadata_is_deep_copied},
      {"owned function data rebinds every pointer",
       test_owned_function_data_rebinds_every_pointer},
      {"leaf tag clone and scope tag copy lifetimes",
       test_leaf_tag_clone_and_scope_tag_copy_lifetimes},
      {"multi-record ancestry uses sibling-local ordinals",
       test_multi_record_ancestry_uses_sibling_local_ordinals},
      {"contradictory multi-record ancestry rejects exact fields",
       test_contradictory_multi_record_ancestry_rejects_exact_fields},
      {"repeated item cannot change from leaf to scope",
       test_repeated_item_cannot_change_from_leaf_to_scope},
      {"new scope item ordinal cannot move backwards across records",
       test_new_scope_item_ordinal_cannot_move_backwards_across_records},
      {"table boundary failures release iterators and clones",
       test_table_boundary_failures_release_iterators_and_clones},
      {"late second target failure discards first target",
       test_late_second_target_failure_discards_first_target},
      {"poison handles are never invoked",
       test_poison_handles_are_never_invoked},
      {"default local and explicit local C GF REAL TL0 pass",
       test_default_local_and_explicit_local_c_gf_real_tl0_pass},
      {"random tag iteration order has equal manifest",
       test_random_tag_iteration_order_has_equal_manifest},
      {"tag flags types counts values and duplicates fail closed",
       test_tag_flags_types_counts_values_and_duplicates_fail_closed},
      {"group and variable full name roundtrip is required",
       test_group_and_variable_full_name_roundtrip_is_required},
      {"scalar array TL1 and unknown region bits reject",
       test_scalar_array_tl1_and_unknown_region_bits_reject},
      {"nonlocal early late singlemap and loop modes reject",
       test_nonlocal_early_late_singlemap_and_loop_modes_reject},
      {"if while storage communication sync and trigger reject",
       test_if_while_storage_communication_sync_and_trigger_reject},
      {"WaveToyX global SYNC is explicitly rejected",
       test_wavetoy_global_sync_is_explicitly_rejected},
      {"raw and normalized access contradiction rejects",
       test_raw_and_normalized_access_contradiction_rejects},
      {"extra missing reordered or changed manifest rejects",
       test_extra_missing_reordered_or_changed_manifest_rejects},
      {"pointer and table handle values do not affect equality",
       test_pointer_and_table_handle_values_do_not_affect_equality},
      {"invalid or mismatched build provenance rejects",
       test_invalid_or_mismatched_build_provenance_rejects},
      {"provenance type sizes match live build and cover observed tags",
       test_provenance_type_sizes_match_live_build_and_cover_observed_tags},
      {"registry publishes only after both targets match",
       test_registry_publishes_only_after_both_targets_match},
      {"native gate observation returns only canonical expectation",
       test_native_gate_observation_returns_only_canonical_expectation},
      {"callback exception is contained and rolls back",
       test_callback_exception_is_contained_and_rolls_back}};
  int failures = 0;
  for (const auto &test : tests)
    failures += run_test(test.first, test.second);
  if (failures != 0)
    std::cerr << failures << " test(s) failed\n";
  else
    std::cout << "Phase 8B3B schedule certification tests passed\n";
  return failures == 0 ? 0 : 1;
}
