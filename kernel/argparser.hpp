#ifndef _ARGPARSER_HPP_
#define _ARGPARSER_HPP_

namespace pic {

  struct CLIArgs {
    std::string journal_file = "journal.txt";
    std::string parameters_file{};
  };

  auto parse_args( int argc, char** argv ) {
    CLIArgs res;

    std::vector<std::string> args;
    {
      args.reserve(argc);
      for ( int i = 0; i < argc; ++i )
        args.push_back({argv[i]});
    }

    for ( int i = 1; i < args.size(); ) {
      if ( "--journal" == args[i] ) {
        res.journal_file = args[i + 1];
        i += 2;
      } else if ( "--params" == args[i] ) {
        res.parameters_file = args[i + 1];
        i += 2;
      }
    }
    return res;
  }
}

#endif
