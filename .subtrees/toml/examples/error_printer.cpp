// This file is a part of toml++ and is subject to the the terms of the MIT license.
// Copyright (c) Mark Gillard <mark.gillard@outlook.com.au>
// See https://github.com/marzer/tomlplusplus/blob/master/LICENSE for the full license text.
// SPDX-License-Identifier: MIT

// This example shows the error messages the library produces by forcing a set of specific parsing
// failures and printing their results.

#include "examples.h"

#define TOML_EXCEPTIONS 0
#define TOML_UNRELEASED_FEATURES 0
#include <toml++/toml.h>

using namespace std::string_view_literals;

namespace
{
	inline constexpr auto invalid_parses = std::array
	{
		"########## comments"sv,
		"# bar\rkek"sv,
		"# bar\bkek"sv,
		"# \xf1\x63"sv,

		"########## inline tables"sv,
		"val = {,}"sv,
		"val = {a='b',}"sv,				// allowed when TOML_UNRELEASED_FEATURES == 1
		"val = {a='b',,}"sv,
		"val = {a='b',"sv,
		"val = {a='b',\n c='d'}"sv,		// allowed when TOML_UNRELEASED_FEATURES == 1
		"val = {?='b'}"sv,

		"########## tables"sv,
		"[]"sv,
		"[foo"sv,
		"[foo] ?"sv,
		"[foo] [bar]"sv,
		"[foo]\n[foo]"sv,
		"? = 'foo' ?"sv,
		"[ [foo] ]"sv,

		"########## arrays"sv,
		"val = [,]"sv,
		"val = ['a',,]"sv,
		"val = ['a',"sv,

		"########## key-value pairs"sv,
		"val = 'foo' ?"sv,
		"val = "sv,
		"val "sv,
		"val ?"sv,
		"val = ]"sv,
		"[foo]\nbar = 'kek'\nbar = 'kek2'"sv,
		"[foo]\nbar = 'kek'\nbar = 7"sv,
		"[foo.bar]\n[foo]\nbar = 'kek'"sv,
		"[foo]\nbar = 'kek'\nbar.kek = 7"sv,
		"[foo]\nbar.? = 'kek'"sv,
		R"('''val''' = 1)"sv,
		R"(a."""val""" = 1)"sv,
		"1= 0x6cA#+\xf1"sv,

		"########## values"sv,
		"val = _"sv,
		"val = G"sv,
		"PATHOLOGICALLY_NESTED"sv, // generated inline

		"########## strings"sv,
		"val = \" \r \""sv,
		R"(val = ")"sv,
		R"(val = "\g")"sv,
		R"(val = "\x20")"sv,			// allowed when TOML_UNRELEASED_FEATURES == 1
		R"(val = "\uFFF")"sv,
		R"(val = "\uFFFG")"sv,
		R"(val = "\UFFFFFFF")"sv,
		R"(val = "\UFFFFFGF")"sv,
		R"(val = "\uD801")"sv,
		R"(val = "\U00110000")"sv,
		R"(val = """ """""")"sv,
		R"(val = ''' '''''')"sv,
		"val = '\n'"sv,

		"########## integers"sv,
		R"(val = -0b0)"sv,
		R"(val = -0o0)"sv,
		R"(val = -0x0)"sv,
		R"(val = +0b0)"sv,
		R"(val = +0o0)"sv,
		R"(val = +0x0)"sv,
		R"(val = 1-)"sv,
		R"(val = -1+)"sv,
		R"(val = -+1)"sv,
		R"(val = 1_0_)"sv,
		R"(val = 1_0_ )"sv,
		R"(val = 999999999999999999999999999999999999 )"sv,
		R"(val = 9223372036854775808 )"sv,
		R"(val = 01 )"sv,

		"########## floats"sv,
		R"(val = 9999999999999999999999999999999999999999999999999999999999999995.0)"sv,
	};

	inline constexpr auto divider =
		"################################################################################"sv;
}

int main()
{
	examples::init();

	for (auto str : invalid_parses)
	{
		if (str.empty())
			continue;

		// section headings
		if (str.substr(0, 10) == "##########"sv)
		{
			std::cout << divider << '\n';
			std::cout << "#    "sv << str.substr(11) << '\n';
			std::cout << divider << "\n\n"sv;
		}

		// error messages
		else
		{
			toml::parse_result result;

			if (str == "PATHOLOGICALLY_NESTED"sv)
			{
				std::string s(1000u, '[');
				constexpr auto start = "array = "sv;
				memcpy(s.data(), start.data(), start.length());
				result = toml::parse(s);
			}
			else
				result = toml::parse(str);

			if (!result)
				std::cout << result.error() << "\n\n"sv;
		}
	}
	return 0;
}
