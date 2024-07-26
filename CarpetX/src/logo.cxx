#include "logo.hxx"

#include <sstream>
#include <string>
#include <vector>

namespace CarpetX {

std::string logo() {

  // buf << "                                                               \n";
  // buf << "    ██████╗ █████╗ ██████╗ ██████╗ ███████╗████████╗██╗  ██╗   \n";
  // buf << "   ██╔════╝██╔══██╗██╔══██╗██╔══██╗██╔════╝╚══██╔══╝╚██╗██╔╝   \n";
  // buf << "   ██║     ███████║██████╔╝██████╔╝█████╗     ██║    ╚███╔╝    \n";
  // buf << "   ██║     ██╔══██║██╔══██╗██╔═══╝ ██╔══╝     ██║    ██╔██╗    \n";
  // buf << "   ╚██████╗██║  ██║██║  ██║██║     ███████╗   ██║   ██╔╝ ██╗   \n";
  // buf << "    ╚═════╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚══════╝   ╚═╝   ╚═╝  ╚═╝   \n";
  // buf << "                                                               \n";

  // colours
  const std::string blk{"\e[30m"};
  const std::string red{"\e[31m"};
  const std::string grn{"\e[32m"};
  const std::string yel{"\e[33m"};
  const std::string blu{"\e[34m"};
  const std::string mag{"\e[35m"};
  const std::string cyn{"\e[36m"};
  const std::string whi{"\e[37m"};
  const std::string reset{"\e[39m"};

  // lines and boxes
  const std::string space{"  "};
  const std::string vline{"▕▏"};
  const std::string hline{"──"};
  const std::string hline2{"━━"};
  const std::string block{"██"};

  const std::vector<std::string> blu_box{
      blu + space + vline + space,
      blu + hline + block + hline,
      blu + space + vline + space,
  };
  const std::vector<std::string> grn_box{
      grn + block + block + block,
      grn + block + block + block,
      grn + block + block + block,
  };
  const std::vector<std::string> red_box{
      red + block + block + block,
      red + block + block + block,
      red + block + block + block,
  };

  const std::vector<std::string> outer_box{
      blu_box[0] + blu_box[0] + blu_box[0],
      blu_box[1] + blu_box[1] + blu_box[1],
      blu_box[2] + blu_box[2] + blu_box[2],
      blu_box[0] + grn_box[0] + blu_box[0],
      blu_box[1] + grn_box[1] + blu_box[1],
      blu_box[2] + grn_box[2] + blu_box[2],
      blu_box[0] + blu_box[0] + blu_box[0],
      blu_box[1] + blu_box[1] + blu_box[1],
      blu_box[2] + blu_box[2] + blu_box[2],
  };
  const std::vector<std::string> inner_box{
      red_box[0] + red_box[0] + red_box[0],
      red_box[1] + red_box[1] + red_box[1],
      red_box[2] + red_box[2] + red_box[2],
      red_box[0] + red_box[0] + red_box[0],
      red_box[1] + red_box[1] + red_box[1],
      red_box[2] + red_box[2] + red_box[2],
      red_box[0] + red_box[0] + red_box[0],
      red_box[1] + red_box[1] + red_box[1],
      red_box[2] + red_box[2] + red_box[2],
  };

  std::vector<std::string> logo;
  for (int line = 0; line < 9; ++line)
    logo.push_back(outer_box[line] + outer_box[line] + outer_box[line]);
  for (int line = 0; line < 9; ++line)
    logo.push_back(outer_box[line] + inner_box[line] + outer_box[line]);
  for (int line = 0; line < 9; ++line)
    logo.push_back(outer_box[line] + outer_box[line] + outer_box[line]);

  std::ostringstream buf;
  buf << "\n";
  for (const auto &str : logo)
    buf << "  " << str << reset << "  \n";
  buf << "\n";
  return buf.str();
}
} // namespace CarpetX
